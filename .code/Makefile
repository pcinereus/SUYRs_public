WS_SRC_DIR  = $(addprefix , .code)
ifdef FILE
WS_FILES = $(WS_SRC_DIR)/$(FILE)
else
WS_FILES ?= $(foreach sdir, $(WS_SRC_DIR), $(wildcard $(sdir)/*.qmd))
endif

HTML_FILES := $(patsubst %.qmd, %.html, $(WS_FILES))

$(info ************************************)
$(info WS Source directory:     $(WS_SRC_DIR))
$(info WS Source files:         $(WS_FILES))
$(info WS HTML files:           $(HTML_FILES))
$(info ************************************)

all: $(HTML_FILES)

$(HTML_FILES) : %.html : %.qmd
	@echo "Compiling documntation for analyses"
	$(info Source = $<; Destination = $@)
	sed -i 's/public\///g' $<
	sed -i 's/,\s*backend\s*=\s*["'\'']cmdstanr?["'\'']\///g' $<
	echo "library(quarto); quarto_render(\"$<\")" | R --no-save --no-restore;

