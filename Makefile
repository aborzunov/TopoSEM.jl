# TopoSEM.jl — reproducibility orchestration.
#
# Common workflows:
#   make deps              install Julia dependencies (both envs)
#   make all               regenerate every data/<sample>_h.csv.gz
#   make <sample>          regenerate one sample (pyramid, sphere, feb_P, ...)
#   make figures           render PNGs in data/ from existing CSVs
#   make visualize         interactive menu — pick a sample and open GLMakie
#   make visualize-<sample> open GLMakie for one sample (e.g. make visualize-pyramid)
#   make test              run the unit-test suite
#   make clean             remove generated CSVs / metadata / PNGs

# `select` builtin lives in bash/zsh; force bash so the menu works regardless
# of the user's default /bin/sh.
SHELL          := /bin/bash

# Two Julia environments:
#   .            — core package + scripts + tests (no OpenGL)
#   examples     — visualisation stack (GLMakie, CairoMakie); split out so
#                  CI / users without a display can avoid pulling ~500 MB of
#                  Makie artefacts.
JULIA          := julia --project=.
JULIA_VIZ      := julia --project=examples
SAMPLES        := jun_vickers jun_vickers_exp1 \
                  jun_sphere jun_sphere_exp2 jun_sphere_exp3 \
                  feb_P feb_PR feb_S feb_V \
                  may_sphere may_piezoceramic may_crooked_60 may_crooked_30
CSV_FILES      := $(addsuffix _h.csv.gz,  $(addprefix data/, $(SAMPLES)))
META_FILES     := $(addsuffix _meta.toml, $(addprefix data/, $(SAMPLES)))
PNG_FILES      := $(addsuffix _topography.png, $(addprefix data/, $(SAMPLES)))
SRC            := $(wildcard src/*.jl) Project.toml Manifest.toml
SCRIPT_COMMON  := scripts/csv_io.jl scripts/reconstruct_common.jl

.PHONY: all clean test figures deps help visualize \
        jun_vickers jun_vickers_exp1 \
        jun_sphere jun_sphere_exp2 jun_sphere_exp3 \
        feb_P feb_PR feb_S feb_V \
        may_sphere may_piezoceramic may_crooked_60 may_crooked_30

# Default: just print help.
help:
	@echo "TopoSEM.jl — reproducibility targets:"
	@echo "  make deps               install Julia dependencies (both envs)"
	@echo "  make all                regenerate all sample reconstructions"
	@echo "  make <sample>           regenerate one of: $(SAMPLES)"
	@echo "  make figures            render PNG figures from existing data/*.csv.gz"
	@echo "  make visualize          interactive menu — pick a sample, open GLMakie"
	@echo "  make visualize-<sample> open GLMakie for one specific sample"
	@echo "  make test               run the unit-test suite"
	@echo "  make clean              remove generated data/*"

all: $(CSV_FILES)

# Pattern rule: every CSV.gz (and its TOML sidecar) depends on the matching
# reconstruction script plus shared infrastructure and the package source.
data/%_h.csv.gz data/%_meta.toml: scripts/reconstruct_%.jl $(SCRIPT_COMMON) $(SRC)
	@mkdir -p data
	$(JULIA) $<

# Override for `jun_sphere_exp3` — documented intentional-failing experiment
# (cubic φ diverges into Inf/NaN around iter ~105, see README §«Polynomial-
# order observation»). The leading `-` makes Make print the failure but NOT
# propagate a non-zero exit, so `make -k all` finishes cleanly after the
# other 12 samples are up-to-date instead of returning Error 2.
data/jun_sphere_exp3_h.csv.gz data/jun_sphere_exp3_meta.toml: \
		scripts/reconstruct_jun_sphere_exp3.jl $(SCRIPT_COMMON) $(SRC)
	@mkdir -p data
	-$(JULIA) $<

# Per-sample shortcuts.
jun_vickers:       data/jun_vickers_h.csv.gz
jun_vickers_exp1:  data/jun_vickers_exp1_h.csv.gz
jun_sphere:        data/jun_sphere_h.csv.gz
jun_sphere_exp2:   data/jun_sphere_exp2_h.csv.gz
jun_sphere_exp3:   data/jun_sphere_exp3_h.csv.gz
feb_P:             data/feb_P_h.csv.gz
feb_PR:            data/feb_PR_h.csv.gz
feb_S:             data/feb_S_h.csv.gz
feb_V:             data/feb_V_h.csv.gz
may_sphere:        data/may_sphere_h.csv.gz
may_piezoceramic:  data/may_piezoceramic_h.csv.gz
may_crooked_60:    data/may_crooked_60_h.csv.gz
may_crooked_30:    data/may_crooked_30_h.csv.gz

# Render PNG figures from the already-shipped CSVs. Implemented as a
# PHONY recipe with an explicit loop (rather than a pattern rule depending
# on `data/%_h.csv.gz`) so that Make never tries to rebuild the CSV if
# its own dependencies (`src/*.jl`, `scripts/reconstruct_common.jl`, …)
# happen to be newer — `make figures` is purely "render existing data".
# Use `make all` or `make <sample>` explicitly when you want fresh CSVs.
figures:
	@for s in $(SAMPLES); do \
	    if [ -f data/$${s}_h.csv.gz ]; then \
	        $(JULIA_VIZ) examples/render_static.jl $$s; \
	    else \
	        echo "skip $$s — data/$${s}_h.csv.gz missing (run 'make $$s' first)"; \
	    fi; \
	done

# Interactive sample picker: shows a numbered menu of all samples, runs
# `examples/visualize.jl <chosen>` over the pre-computed CSV in `data/`.
# `select` is a bash builtin, hence `SHELL := /bin/bash` at the top.
visualize:
	@PS3=$$'\nPick a sample number (q to quit): '; \
	echo "Available samples:"; \
	select s in $(SAMPLES); do \
	    case "$$REPLY" in \
	        q|Q) echo "cancelled."; break ;; \
	        *) if [ -n "$$s" ]; then \
	               $(JULIA_VIZ) examples/visualize.jl "$$s"; \
	               break; \
	           else \
	               echo "invalid choice — try again or 'q' to quit"; \
	           fi ;; \
	    esac; \
	done

# Direct shortcut: `make visualize-pyramid` opens GLMakie for `pyramid` straight away.
visualize-%:
	$(JULIA_VIZ) examples/visualize.jl $*

test:
	$(JULIA) -e 'using Pkg; Pkg.test()'

deps:
	$(JULIA)     -e 'using Pkg; Pkg.instantiate()'
	$(JULIA_VIZ) -e 'using Pkg; Pkg.instantiate()'

clean:
	rm -f $(CSV_FILES) $(META_FILES) $(PNG_FILES)
