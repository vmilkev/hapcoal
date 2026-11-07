#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <cassert>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum class ModelType { SMC, SMC_Prime, MaCS };

struct ChromosomeConfig {
    int id = 1;
    double length = 200000.0;       // Physical size bounds (200 kb scalable to hundreds of Mb)
    double mu = 1.5e-7;             // Calibrated mutation density to prevent small-footprint variant noise
    double rho_background = 1e-8;   // Baseline physical crossover rate per bp per generation
    std::string genetic_map_path = "";
    double centromere_start_pct = 0.45; 
    double centromere_end_pct = 0.55;   
    double default_hotspot_density_per_mb = 20.0; 
};

struct Options {
    int num_samples = 100;             // Fully modular for thousands of production samples
    double N0 = 10000.0;               // Scaling population parameters
    ModelType model = ModelType::MaCS;  
    double macs_history_cutoff = 3.5;  // Bounded Markov history window
    std::vector<ChromosomeConfig> chromosomes;
};

struct Hotspot {
    double center;
    double amplitude_multiplier; 
    double width;
};

struct MapEntry {
    int chr_idx;         
    double physical_pos; 
    double rate;         
    double genetic_pos;  
};

struct CoalescentNode {
    int id = -1;
    int parent = -1;
    int left = -1;
    int right = -1;
    double time = 0.0;
};

struct TransposedSNP {
    double physical_pos;
    size_t matrix_offset;
};

class LocalTree {
public:
    std::vector<CoalescentNode> nodes;
    int root_id = -1;
    int num_samples = 0;
    int next_node_id = 0;
    double total_length = 0.0; 
    std::vector<int> active_branches;     
    std::vector<int> node_to_branch_idx;  

    LocalTree(int n) : num_samples(n) {
        size_t max_nodes = static_cast<size_t>(8 * n + 10000); 
        nodes.resize(max_nodes); 
        node_to_branch_idx.resize(max_nodes, -1);
        clear_pool();
    }

    void clear_pool() {
        next_node_id = num_samples;
        root_id = -1;
        total_length = 0.0;
        active_branches.clear();
        std::fill(node_to_branch_idx.begin(), node_to_branch_idx.end(), -1);

        for (size_t i = 0; i < nodes.size(); ++i) {
            nodes[i].id = static_cast<int>(i);
            nodes[i].parent = -1;
            nodes[i].left = -1;
            nodes[i].right = -1;
            nodes[i].time = 0.0;
        }
    }

    void register_branch(int node_id) {
        if (node_id == root_id || nodes[node_id].parent == -1) {
            unregister_branch(node_id);
            return;
        }
        if (node_to_branch_idx[node_id] == -1) {
            node_to_branch_idx[node_id] = static_cast<int>(active_branches.size());
            active_branches.push_back(node_id);
            
            double branch_len = nodes[nodes[node_id].parent].time - nodes[node_id].time;
            if (branch_len > 0.0) total_length += branch_len;
        }
    }

    void unregister_branch(int node_id) {
        int idx = node_to_branch_idx[node_id];
        if (idx == -1) return;

        double branch_len = nodes[nodes[node_id].parent].time - nodes[node_id].time;
        if (branch_len > 0.0) total_length -= branch_len;

        int last_node_id = active_branches.back();
        active_branches[idx] = last_node_id;
        node_to_branch_idx[last_node_id] = idx;

        active_branches.pop_back();
        node_to_branch_idx[node_id] = -1;
    }
    void build_initial_kingman(std::mt19937_64& rng) {
        clear_pool();
        std::vector<int> active_lineages;
        active_lineages.reserve(num_samples);
        for (int i = 0; i < num_samples; ++i) {
            active_lineages.push_back(i);
        }

        double current_time = 0.0;
        while (active_lineages.size() > 1) {
            int k = static_cast<int>(active_lineages.size());
            double rate = k * (k - 1) / 2.0;
            std::exponential_distribution<double> exp_dist(rate);
            current_time += exp_dist(rng);

            std::uniform_int_distribution<int> uni_dist(0, k - 1);
            int idx1 = uni_dist(rng);
            int idx2 = uni_dist(rng);
            while (idx1 == idx2) idx2 = uni_dist(rng);

            int u = active_lineages[idx1];
            int v = active_lineages[idx2];

            int parent_idx = next_node_id++;
            assert(parent_idx < static_cast<int>(nodes.size()) && "LocalTree node pool limit reached");
            
            nodes[parent_idx].time = current_time;
            nodes[parent_idx].left = u;
            nodes[parent_idx].right = v;
            nodes[u].parent = parent_idx;
            nodes[v].parent = parent_idx;

            if (idx1 < idx2) std::swap(idx1, idx2);
            
            active_lineages[idx1] = active_lineages.back();
            active_lineages.pop_back();
            
            active_lineages[idx2] = active_lineages.back();
            active_lineages.pop_back();
            
            active_lineages.push_back(parent_idx);
        }
        root_id = active_lineages[0];

        total_length = 0.0;
        for (int i = 0; i < next_node_id; ++i) {
            register_branch(i);
        }
    }
};

double genetic_to_physical(double gen_pos, const std::vector<MapEntry>& active_map, double chr_length) {
    if (active_map.empty() || gen_pos <= active_map.front().genetic_pos) return 0.0;
    if (gen_pos >= active_map.back().genetic_pos) return chr_length;

    auto it = std::upper_bound(active_map.begin(), active_map.end(), gen_pos, [](double val, const MapEntry& entry) {
        return val < entry.genetic_pos;
    });

    auto prev = std::prev(it);
    double g_diff = it->genetic_pos - prev->genetic_pos;
    if (g_diff <= 0.0) return prev->physical_pos;

    double fraction = (gen_pos - prev->genetic_pos) / g_diff;
    return prev->physical_pos + fraction * (it->physical_pos - prev->physical_pos);
}

double physical_to_genetic(double phys_pos, const std::vector<MapEntry>& active_map) {
    if (active_map.empty() || phys_pos <= active_map.front().physical_pos) return 0.0;
    if (phys_pos >= active_map.back().physical_pos) return active_map.back().genetic_pos;

    auto it = std::upper_bound(active_map.begin(), active_map.end(), phys_pos, [](double val, const MapEntry& entry) {
        return val < entry.physical_pos;
    });

    auto prev = std::prev(it);
    double p_diff = it->physical_pos - prev->physical_pos;
    if (p_diff <= 0.0) return prev->genetic_pos;

    double fraction = (phys_pos - prev->physical_pos) / p_diff;
    return prev->genetic_pos + fraction * (it->genetic_pos - prev->genetic_pos);
}

void collect_mutated_leaves_iterative(int start_node_id, const LocalTree& tree, std::vector<int>& leaves) {
    std::vector<int> stack_buffer; 
    stack_buffer.clear();
    stack_buffer.push_back(start_node_id);

    while (!stack_buffer.empty()) {
        int curr = stack_buffer.back();
        stack_buffer.pop_back();

        if (curr < tree.num_samples) {
            leaves.push_back(curr);
        } else {
            if (tree.nodes[curr].right != -1) stack_buffer.push_back(tree.nodes[curr].right);
            if (tree.nodes[curr].left != -1) stack_buffer.push_back(tree.nodes[curr].left);
        }
    }
}

void find_active_lineages_iterative(int start_node_id, double time, const LocalTree& tree, std::vector<int>& lineages) {
    if (start_node_id == -1) return;
    std::vector<int> stack_buffer;
    stack_buffer.clear();
    stack_buffer.push_back(start_node_id);

    while (!stack_buffer.empty()) {
        int curr = stack_buffer.back();
        stack_buffer.pop_back();

        if (tree.nodes[curr].time <= time) {
            if (tree.nodes[curr].parent == -1 || tree.nodes[tree.nodes[curr].parent].time > time) {
                lineages.push_back(curr);
                continue;
            }
        }
        if (tree.nodes[curr].right != -1) stack_buffer.push_back(tree.nodes[curr].right);
        if (tree.nodes[curr].left != -1) stack_buffer.push_back(tree.nodes[curr].left);
    }
}
std::vector<MapEntry> load_genetic_map(const std::string& path) {
    return std::vector<MapEntry>(); // Placeholder stub for compilation context
}

std::vector<MapEntry> generate_hotspot_map(const ChromosomeConfig& chr, const std::vector<Hotspot>& hotspots) {
    std::vector<MapEntry> map_entries;
    int steps = 5000; 
    double step_size = chr.length / steps;
    double accumulated_genetic = 0.0;

    map_entries.push_back({chr.id, 0.0, chr.rho_background, 0.0});

    for (int i = 1; i <= steps; ++i) {
        double phys_pos = i * step_size;
        double current_rate = chr.rho_background;

        double multiplier = 1.0;
        for (const auto& hs : hotspots) {
            double dist = phys_pos - hs.center;
            multiplier += (hs.amplitude_multiplier - 1.0) * std::exp(-(dist * dist) / (2.0 * hs.width * hs.width));
        }
        current_rate *= multiplier;

        accumulated_genetic += current_rate * step_size;
        map_entries.push_back({chr.id, phys_pos, current_rate, accumulated_genetic});
    }
    return map_entries;
}

int select_branch_by_cdf(const LocalTree& tree, std::mt19937_64& rng) {
    if (tree.active_branches.empty()) return tree.root_id;

    std::uniform_real_distribution<double> dist(0.0, tree.total_length);
    double target = dist(rng);
    
    double accumulated = 0.0;
    for (int branch_id : tree.active_branches) {
        double branch_len = tree.nodes[tree.nodes[branch_id].parent].time - tree.nodes[branch_id].time;
        if (branch_len > 0.0) {
            accumulated += branch_len;
            if (target <= accumulated) {
                return branch_id;
            }
        }
    }
    return tree.active_branches.back();
}

void simulate_and_stream_chromosome(
    const Options& opts, 
    const ChromosomeConfig& chr, 
    const std::vector<MapEntry>& active_map, 
    std::mt19937_64& rng,
    std::ofstream& hap_file,
    std::ofstream& ld_file
) {
    LocalTree tree(opts.num_samples);
    tree.build_initial_kingman(rng);
    
    std::vector<int> free_node_ids;
    free_node_ids.reserve(tree.nodes.size());
    for (size_t i = tree.next_node_id; i < tree.nodes.size(); ++i) {
        free_node_ids.push_back(static_cast<int>(i));
    }

    double current_phys_pos = 0.0;
    std::vector<TransposedSNP> snp_metadata; 
    
    // OPTIMIZATION 2: Pre-calculate expected mutations and pre-allocate exactly
    // Expected mutations ~ chr.length * chr.mu * 2*N0 * 4*log(n) / 2
    double expected_branch_len = 4.0 * std::log(opts.num_samples);
    size_t expected_snps = static_cast<size_t>(chr.length * chr.mu * 2.0 * opts.N0 * expected_branch_len) + 1000;
    expected_snps = std::min(expected_snps, size_t(chr.length / 10)); // cap at 1 SNP per 10bp
    size_t words_per_sample = (opts.num_samples + 63) / 64;
    std::vector<uint64_t> flat_haplotype_matrix;
    flat_haplotype_matrix.reserve(expected_snps * words_per_sample);
    snp_metadata.reserve(expected_snps);

    std::vector<int> leaf_buffer;
    leaf_buffer.reserve(opts.num_samples);

    // --- MAIN SPATIAL PROGRESSION WINDOW LOOP ---
    while (current_phys_pos < chr.length) {
        double total_branch_length = tree.total_length;
        if (total_branch_length <= 0.0) break;

        double current_rho = chr.rho_background;
        auto it = std::upper_bound(active_map.begin(), active_map.end(), current_phys_pos, 
            [](double val, const MapEntry& entry) { return val < entry.physical_pos; });
        if (it != active_map.end() && it != active_map.begin()) {
            current_rho = std::prev(it)->rate;
        }

        // Biological spatial scaling step generated directly in physical base pairs (bp)
        double physical_recombination_rate = total_branch_length * (2.0 * opts.N0) * current_rho;
        double delta_x = chr.length - current_phys_pos;
        
        if (physical_recombination_rate > 0.0) {
            std::exponential_distribution<double> exp_dist(physical_recombination_rate);
            double sampled_delta = exp_dist(rng);
            if (sampled_delta < delta_x) {
                delta_x = sampled_delta;
            }
        }

        if (delta_x < 1.0) delta_x = 1.0;
        double next_phys_pos = current_phys_pos + delta_x;
        if (next_phys_pos > chr.length) next_phys_pos = chr.length;
        
        double physical_span = next_phys_pos - current_phys_pos;
        // --- SPATIAL MUTATION GENERATION WINDOW ---
        if (physical_span > 0.0 && !tree.active_branches.empty()) {
            double total_branch_length_generations = total_branch_length * (2.0 * opts.N0);
            double total_chunk_expected_mu = chr.mu * total_branch_length_generations * physical_span;
            
            std::poisson_distribution<int> poisson_muts(total_chunk_expected_mu);
            int muts_in_chunk = poisson_muts(rng);

            if (muts_in_chunk > 0) {
                std::uniform_real_distribution<double> pos_in_chunk(current_phys_pos, next_phys_pos);

                for (int m = 0; m < muts_in_chunk; ++m) {
                    double snp_pos = pos_in_chunk(rng);
                    
                    // FIX: Select target mutation branch proportional to its length using the CDF method
                    int target_mut_node = select_branch_by_cdf(tree, rng);

                    size_t offset = flat_haplotype_matrix.size();
                    flat_haplotype_matrix.resize(offset + words_per_sample, 0ULL);

                    leaf_buffer.clear();
                    collect_mutated_leaves_iterative(target_mut_node, tree, leaf_buffer);
                    for (int leaf : leaf_buffer) {
                        flat_haplotype_matrix[offset + (leaf / 64)] |= (1ULL << (leaf % 64));
                    }
                    snp_metadata.push_back({snp_pos, offset});
                }
            }
        }

        if (next_phys_pos >= chr.length) {
            current_phys_pos = chr.length;
            break;
        }

        // --- TREE RECOMBINATION DISSECTING SELECTION ---
        int broken_node = select_branch_by_cdf(tree, rng);

        double low_time = tree.nodes[broken_node].time;
        double high_time = tree.nodes[tree.nodes[broken_node].parent].time;
        std::uniform_real_distribution<double> height_dist(low_time, high_time);
        double recomb_time = height_dist(rng);

        int old_parent = tree.nodes[broken_node].parent;
        tree.unregister_branch(broken_node);
        tree.unregister_branch(tree.nodes[old_parent].left);
        tree.unregister_branch(tree.nodes[old_parent].right);

        if (tree.nodes[old_parent].left == broken_node) tree.nodes[old_parent].left = -1;
        else tree.nodes[old_parent].right = -1;

        int sister_node = (tree.nodes[old_parent].left == -1) ? tree.nodes[old_parent].right : tree.nodes[old_parent].left;
        int grandparent = tree.nodes[old_parent].parent;
        tree.unregister_branch(sister_node);

        if (grandparent != -1) {
            tree.unregister_branch(old_parent);
            if (tree.nodes[grandparent].left == old_parent) tree.nodes[grandparent].left = sister_node;
            else tree.nodes[grandparent].right = sister_node;
            tree.nodes[sister_node].parent = grandparent;
            tree.register_branch(sister_node);
            if (tree.nodes[grandparent].parent != -1) tree.register_branch(grandparent);
        } else {
            tree.nodes[sister_node].parent = -1; 
            tree.root_id = sister_node;
        }

        tree.nodes[old_parent].parent = -1; 
        tree.nodes[old_parent].left = -1; 
        tree.nodes[old_parent].right = -1; 
        tree.nodes[old_parent].time = 0.0;
        free_node_ids.push_back(old_parent);

        // --- FIXED ACTIVE LINEAGE EXTRACTION (DUAL-FRAGMENT TRACKING) ---
        std::vector<int> active_lineages;
        find_active_lineages_iterative(tree.root_id, recomb_time, tree, active_lineages);
        find_active_lineages_iterative(broken_node, recomb_time, tree, active_lineages);
        
        // --- MULTI-MODEL STRUCTURAL COALESCENCE FILTERS ---
        if (opts.model == ModelType::SMC) {
            active_lineages.erase(
                std::remove_if(active_lineages.begin(), active_lineages.end(), 
                    [&](int lineage_id) {
                        return (lineage_id == sister_node && recomb_time < tree.nodes[old_parent].time);
                    }), 
                active_lineages.end()
            );
        } 
        else if (opts.model == ModelType::MaCS) {
            active_lineages.erase(
                std::remove_if(active_lineages.begin(), active_lineages.end(), 
                    [&](int lineage_id) {
                        return (tree.nodes[lineage_id].time > (recomb_time + opts.macs_history_cutoff));
                    }), 
                active_lineages.end()
            );
        }
        else if (opts.model == ModelType::SMC_Prime) {
            active_lineages.erase(
                std::remove_if(active_lineages.begin(), active_lineages.end(),
                    [&](int lineage_id) {
                        return (lineage_id == sister_node && recomb_time < tree.nodes[old_parent].time);
                    }), 
                active_lineages.end()
            );
        }

        // --- LINEAGE GRAFTING MECHANICS ---
        double coal_time = recomb_time;
        int target_coal_node = -1;
        if (active_lineages.empty()) {
            std::exponential_distribution<double> top_dist(1.0);
            double calculated_top = tree.nodes[tree.root_id].time + top_dist(rng);
            if (opts.model == ModelType::MaCS) {
                double max_allowable = recomb_time + opts.macs_history_cutoff;
                coal_time = (calculated_top > max_allowable) ? max_allowable : calculated_top;
            } else {
                coal_time = calculated_top;
            }
            target_coal_node = tree.root_id;
        } else {
            int k = static_cast<int>(active_lineages.size());
            std::exponential_distribution<double> exp_coal(k);
            coal_time += exp_coal(rng);
            
            if (opts.model == ModelType::MaCS && coal_time > (recomb_time + opts.macs_history_cutoff)) {
                coal_time = recomb_time + opts.macs_history_cutoff;
            }

            std::uniform_int_distribution<int> pick_dist(0, k - 1);
            target_coal_node = active_lineages[pick_dist(rng)];
            while (tree.nodes[target_coal_node].parent != -1 && coal_time > tree.nodes[tree.nodes[target_coal_node].parent].time) {
                target_coal_node = tree.nodes[target_coal_node].parent;
            }
        }

        tree.unregister_branch(target_coal_node);
        int new_parent = free_node_ids.back(); 
        free_node_ids.pop_back();

        tree.nodes[new_parent].time = coal_time; 
        tree.nodes[new_parent].left = broken_node; 
        tree.nodes[new_parent].right = target_coal_node;
        
        int target_parent = tree.nodes[target_coal_node].parent;
        tree.nodes[broken_node].parent = new_parent; 
        tree.nodes[target_coal_node].parent = new_parent; 
        tree.nodes[new_parent].parent = target_parent;

        if (target_parent != -1) {
            tree.unregister_branch(target_parent);
            if (tree.nodes[target_parent].left == target_coal_node) tree.nodes[target_parent].left = new_parent;
            else tree.nodes[target_parent].right = new_parent;
        } else {
            tree.root_id = new_parent;
        }

        tree.register_branch(broken_node); 
        tree.register_branch(target_coal_node);
        tree.register_branch(new_parent);
        if (target_parent != -1) tree.register_branch(target_parent);
        current_phys_pos = next_phys_pos;
    }
    // --- SORTING AND FILE STREAM DATA OUTPUTS SERIALIZATION ---
    std::sort(snp_metadata.begin(), snp_metadata.end(),
        [](const TransposedSNP& a, const TransposedSNP& b) { return a.physical_pos < b.physical_pos; });
    
    int total_snps = static_cast<int>(snp_metadata.size());
    std::cout << " ↳ Generated " << total_snps << " valid active SNPs for Chromosome [" << chr.id << "]." << std::endl;
    
    for (int i = 0; i < opts.num_samples; ++i) {
        hap_file << chr.id << ",Sample_" << i << ",";
        for (int s = 0; s < total_snps; ++s) {
            size_t offset = snp_metadata[s].matrix_offset;
            bool bit_val = (flat_haplotype_matrix[offset + (i / 64)] & (1ULL << (i % 64))) != 0ULL;
            hap_file << (bit_val ? '1' : '0');
        }
        hap_file << "\n";
    }

    std::vector<double> p_A(total_snps, 0.0);
    for (int s = 0; s < total_snps; ++s) {
        int counts = 0;
        size_t offset = snp_metadata[s].matrix_offset;
        for (size_t w = 0; w < words_per_sample; ++w) {
            uint64_t v = flat_haplotype_matrix[offset + w];
            for (; v; ++counts) { v &= v - 1; }
        }
        p_A[s] = static_cast<double>(counts) / opts.num_samples;
    }

// --- HIGH-SPEED PRODUCTION SLIDING WINDOW LD SCANNER ---
    const double bin_size = 1000.0;
    const double max_ld_distance = 200000.0;

    // OPTIMIZATION 1: Binary search on sorted SNP positions
    for (int i = 0; i < total_snps; ++i) {
        size_t offset_i = snp_metadata[i].matrix_offset;
        if (p_A[i] <= 0.05 || p_A[i] >= 0.95) continue;

        // Binary search to find the last SNP within max_ld_distance
        int left_j = i + 1;
        int right_j = total_snps - 1;
        int upper_j = i;
        while (left_j <= right_j) {
            int mid = left_j + (right_j - left_j) / 2;
            double d = snp_metadata[mid].physical_pos - snp_metadata[i].physical_pos;
            if (d > max_ld_distance) {
                right_j = mid - 1;
            } else {
                upper_j = mid;
                left_j = mid + 1;
            }
        }

        for (int j = i + 1; j <= upper_j; ++j) {
            double distance = snp_metadata[j].physical_pos - snp_metadata[i].physical_pos;
            if (p_A[j] <= 0.05 || p_A[j] >= 0.95) continue;

            int shared_counts = 0;
            size_t offset_j = snp_metadata[j].matrix_offset;

            for (size_t w = 0; w < words_per_sample; ++w) {
                uint64_t intersection = flat_haplotype_matrix[offset_i + w] & flat_haplotype_matrix[offset_j + w];
                for (; intersection; ++shared_counts) { intersection &= intersection - 1; }
            }

            double p_AB = static_cast<double>(shared_counts) / opts.num_samples;
            double D = p_AB - (p_A[i] * p_A[j]);
            double denom = p_A[i] * (1.0 - p_A[i]) * p_A[j] * (1.0 - p_A[j]);
            if (denom <= 0.0) continue;

            double r2 = (D * D) / denom;
            int current_bin = static_cast<int>(std::floor(distance / bin_size)) * static_cast<int>(bin_size);
            ld_bins[current_bin].first += r2;
            ld_bins[current_bin].second += 1;
        }
    }

    for (int b = 0; b < num_bins; ++b) {
        if (ld_bins_obs[b] > 0) {
            ld_file << chr.id << "," << (b * bin_size) << "," << std::fixed << std::setprecision(5) << (ld_bins_sum[b] / ld_bins_obs[b]) << "\n";
        }
    }
}

int main() {
    std::mt19937_64 rng(42); 
    Options opts;
    opts.num_samples = 1000; // Ready for production scaling up to thousands of samples
    opts.N0 = 10000.0;
    
    opts.model = ModelType::MaCS; 
    //opts.model = ModelType::SMC; 
    //opts.model = ModelType::SMC_Prime; 
    opts.macs_history_cutoff = 3.5; 

    // Production-ready Baseline Calibration setup:
    double length_chr1 = 200000000.0; 
    double length_chr2 = 200000000.0; 
    
    // Baseline physical recombination rate per base pair per generation (e.g. 1 cM/Mb)
    double true_rho_background = 1.0e-8; 

    ChromosomeConfig chr1{1, length_chr1, 1.5e-8, true_rho_background, "", 0.45, 0.55, 20.0};
    ChromosomeConfig chr2{2, length_chr2, 1.5e-8, true_rho_background, "", 0.40, 0.50, 20.0};
    opts.chromosomes.push_back(chr1);
    opts.chromosomes.push_back(chr2);

    std::ofstream hap_file("haplotypes.csv");
    hap_file << "Chromosome,Sample_ID,Sequence_Haplotype\n";

    std::ofstream ld_file("ld_decay.csv");
    ld_file << "Chromosome,Distance_Bin_bp,Mean_R2_Linkage_Disequilibrium\n";

    std::ofstream calculated_map_file("calculated_linkage_map.csv");
    calculated_map_file << "Chromosome,Physical_Position_bp,Recombination_Rate,Estimated_Genetic_Position_cM\n";

    std::cout << "🧬 Starting Multi-Chromosome Genome Production System (Optimized & Calibrated)..." << std::endl;

    for (const auto& chr : opts.chromosomes) {
        std::cout << "\nParsing Chromosome [" << chr.id << "] | Bounds: " << chr.length << " bp" << std::endl;
        
        std::vector<MapEntry> active_map;
        bool user_submitted_map = (!chr.genetic_map_path.empty());

        if (user_submitted_map) {
            active_map = load_genetic_map(chr.genetic_map_path);
        }

        if (active_map.empty()) {
            std::vector<Hotspot> sampled_hotspots;
            double chr_length_mb = chr.length / 1e6;
            int total_hotspots_to_sample = static_cast<int>(std::round(chr.default_hotspot_density_per_mb * chr_length_mb));
            if (total_hotspots_to_sample < 2) total_hotspots_to_sample = 2;

            double centromere_start_bp = chr.length * chr.centromere_start_pct;
            double centromere_end_bp = chr.length * chr.centromere_end_pct;

            std::uniform_real_distribution<double> position_dist(0.0, chr.length);
            // Hotspots act as standard 10x-20x regional multipliers over background recombination levels
            std::uniform_real_distribution<double> amplitude_dist(10.0, 20.0); 
            std::uniform_real_distribution<double> width_dist(800.0, 1500.0);

            while (static_cast<int>(sampled_hotspots.size()) < total_hotspots_to_sample) {
                double candidate_center = position_dist(rng);
                if (candidate_center >= centromere_start_bp && candidate_center <= centromere_end_bp) continue;
                sampled_hotspots.push_back({candidate_center, amplitude_dist(rng), width_dist(rng)});
            }
            active_map = generate_hotspot_map(chr, sampled_hotspots);
        }

        simulate_and_stream_chromosome(opts, chr, active_map, rng, hap_file, ld_file);

        if (!user_submitted_map) {
            for (const auto& entry : active_map) {
                calculated_map_file << chr.id << "," << std::fixed << std::setprecision(2) << entry.physical_pos << ","
                                    << std::setprecision(5) << entry.rate << "," << std::setprecision(4) << (entry.genetic_pos * 100.0) << "\n";
            }
        }
    }

    hap_file.close(); ld_file.close(); calculated_map_file.close();
    std::cout << "\n✅ Execution successful! Calibrated decay profiles written directly to disk.\n";
    return 0;
}
