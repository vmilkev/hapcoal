/**
 * @file hapcoal.hpp
 * @brief Umbrella header of the hapcoal library — include this from client code.
 *
 * Link with `lib/libhapcoal.a` plus the OpenMP runtime and zlib
 * (`make print-flags` prints the exact compiler and linker flags).
 *
 * Public components:
 *  - CoalescentSimulator.hpp : parameters, results, SMC / SMC′ / MaCS simulator
 *  - Demography.hpp          : piecewise-constant population size history
 *  - GeneticMap.hpp          : genetic maps (files, hotspot landscapes)
 *  - Haplotypes.hpp          : packed haplotype containers, `.hap.bin` reader
 *  - LDCalculator.hpp        : r², LD decay, LD matrices
 *  - MapEstimator.hpp        : realised and LD-based genetic maps
 *  - VcfReader.hpp           : phased VCF input, observed LD
 *  - DemographyFit.hpp       : fit a demography to an observed LD decay
 */
#pragma once

#define HAPCOAL_VERSION_MAJOR 1
#define HAPCOAL_VERSION_MINOR 0
#define HAPCOAL_VERSION_PATCH 0
#define HAPCOAL_VERSION_STRING "1.0.0"

#include "hapcoal/CoalescentSimulator.hpp"
#include "hapcoal/Demography.hpp"
#include "hapcoal/DemographyFit.hpp"
#include "hapcoal/GeneticMap.hpp"
#include "hapcoal/Haplotypes.hpp"
#include "hapcoal/LDCalculator.hpp"
#include "hapcoal/MapEstimator.hpp"
#include "hapcoal/VcfReader.hpp"
