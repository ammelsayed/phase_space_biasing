include ../../make_opts

all: ptl1_met_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module ptl1_met_bias
#
ptl1_met_bias.o: ptl1_met_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o ptl1_met_bias.o ptl1_met_bias.f

ptl1_met_bias: ptl1_met_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:
	@echo "VALID"
