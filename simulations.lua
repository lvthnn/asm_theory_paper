-------------------------------------------------------------------------------
-- SIMULATIONS FOR 'A HIGH-DIMENSIONAL PERSPECTIVE OF ASSORTATIVE MATING'
-------------------------------------------------------------------------------
-- To run these simulations, you will need a working installation of knoop. See
-- https://github.com/knoop-dev/knoop.
--
-- Some across-scenario simulation parameters. You may want to tweak the
-- multithreading options depending whether you have less/more available threads
-- on your machine, and/or the number of replicate simulations performed, etc.
-- More replicates certainly never hurt anyone (excepting the events of Blade
-- Runner).

local params = {
	n_replicates = 50,
	n_workers = 3,
	n_threads_per_worker = 4,
	n_individuals = 100000,
	n_loci = 50000,
	n_generations = 50,
	locus_freq = beta { shape1 = 4, shape2 = 4 },
	locus_rec = 0.50,
	locus_mut = 0.00,
	output_dir = "results",
}

local random = uniform {}

--- Generate a random orthogonal matrix
---
---@param dim integer Number of rows / columns in the matrix
---
---@return matrix
local function generate_orthogonal(dim)
	local U, _ = qr(matrix { random(dim * dim), dim, dim })
	return U
end

--- Generate generatic traits for the simulation
---
---@param n_traits integer Number of traits to generate
---@param var_g matrix Genetic component covariance
---@param var_e matrix Environmental component covariance
---
---@return table
local function generate_traits(n_traits, var_g, var_e)
	local traits = {}

	for j = 1, n_traits do
		traits[j] = phenotype {
			name = "trait" .. j,
			n_loci = params.n_loci // n_traits,
			var_genetic = var_g[j],
			var_environmental = var_e[j],
			effects = normal {},
		}
	end

	return traits
end

--- Run a knoop simulation
---
---@param spec_g vector Genetic component covariance spectrum
---@param spec_e vector Environmental component covariance spectrum
---@param spec_m vector (Symmetrised) assorting matrix spectrum
---@param commuting boolean Share eigenbasis for the above?
---@param rng_seed integer Random seed used for this simulation
---@param name string Descriptive name of this simulation config
local function run_scenario(spec_g, spec_e, spec_m, commuting, rng_seed, name)
	print("Running scenario " .. name)

	commuting = commuting or false

	-- Sternly reprimand the user if they misbehave
	local n_traits = #spec_g

	if #spec_e ~= n_traits then
		error "Spectrum of spec_e must be of same length as spec_g"
	end

	if #spec_m ~= n_traits then
		error "Spectrum of spec_m must be of same length as spec_g"
	end

	seed(rng_seed)

	-- If commuting is enabled, we use the same eigenbasis (i.e., the same matrix
	-- U) for all of the matrices cov_g, cov_e, and mate_cor.
	local U_g, U_e, U_m

	if commuting then
		local U = generate_orthogonal(n_traits)
		U_g = U
		U_e = U
		U_m = U
	else
		U_g = generate_orthogonal(n_traits)
		U_e = generate_orthogonal(n_traits)
		U_m = generate_orthogonal(n_traits)
	end

	-- Set up the covariance matrices and assorting matrix. Since knoop takes the
	-- panmictic genetic and environmental component covariances as correlation
	-- matrices, must multiply the inverse standard deviations.
	local cov_g = U_g * diag(spec_g) * U_g.T
	local cov_e = U_e * diag(spec_e) * U_e.T
	local cor_m = U_m * diag(spec_m) * U_m.T

	local diag_inv_g = diag(1 / sqrt(diag(cov_g)))
	local diag_inv_e = diag(1 / sqrt(diag(cov_e)))

	local cor_g = diag_inv_g * cov_g * diag_inv_g
	local cor_e = diag_inv_e * cov_e * diag_inv_e

	-- Alrighty. We are now ready to cobble together the simulation.
	local genome = genome {
		n_loci = params.n_loci,
		locus_freq = params.locus_freq,
		locus_rec = params.locus_rec,
		locus_mut = params.locus_mut,
	}

	local phenotypes = phenotypes {
		phenotypes = generate_traits(n_traits, diag(cov_g), diag(cov_e)),
		gen_cor = cor_g,
		env_cor = cor_e,
	}

	local mating = mating {
		type = mating_type.assortative,
		cor = cor_m,
		whiten = true,
		tol = 0.025,
	}

	local estimators = {
		estimator.heritability {},
		estimator.phenotype_cov { component.genetic },
		estimator.phenotype_cov { component.total },
		estimator.spouse_cor { component.total },
		estimator.sibling_cov { component.genetic },
		estimator.cousin_cov { 1, component.genetic },
		estimator.cousin_cov { 2, component.genetic },
		estimator.cousin_cov { 3, component.genetic },
		estimator.ancestor_cov { 1, component.genetic },
		estimator.ancestor_cov { 2, component.genetic },
	}

	simulation {
		n_individuals = params.n_individuals,
		n_generations = params.n_generations,
		n_replicates = params.n_replicates,
		n_workers = params.n_workers,
		n_threads_per_worker = params.n_threads_per_worker,
		output_dir = params.output_dir,
		output_name = name,
		genome = genome,
		phenotypes = phenotypes,
		mating = mating,
		estimators = estimators,
		pedigree_max_depth = 5,
		pedigree_warmup = true,
		log_to_file = false,
		log_level = log_level.info,
	}
end

-- RUN THE SIMULATIONS

local negative = vector { -1, -1, -1, -1, -1 }
local indefinite = vector { 1, -1, 1, -1, 1 }
local spec_g = vector { 0.50, 0.40, 0.35, 0.30, 0.25 }
local spec_e = vector { 0.50, 0.60, 0.65, 0.70, 0.75 }
local spec_m = vector { 0.60, 0.45, 0.30, 0.15, 0.05 }

run_scenario(spec_g, spec_e, spec_m, true, 25783571, "asm_positive_commuting")

run_scenario(
	spec_g,
	spec_e,
	spec_m,
	false,
	18219857,
	"asm_positive_noncommuting"
)

run_scenario(
	spec_g,
	spec_e,
	negative * spec_m,
	true,
	958109218,
	"asm_negative_commuting"
)

run_scenario(
	spec_g,
	spec_e,
	negative * spec_m,
	false,
	75117114,
	"asm_negative_noncommuting"
)

run_scenario(
	spec_g,
	spec_e,
	indefinite * spec_m,
	true,
	5971978,
	"asm_indefinite_commuting"
)

run_scenario(
	spec_g,
	spec_e,
	indefinite * spec_m,
	false,
	48588215,
	"asm_indefinite_noncommuting"
)
