include ../../make_opts

all: met_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module met_bias
#
met_bias.o: met_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o met_bias.o met_bias.f

met_bias: met_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:
	@echo "VALID"
