include ../../make_opts

all: ptl1_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module ptl1_bias
#
ptl1_bias.o: ptl1_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o ptl1_bias.o ptl1_bias.f

ptl1_bias: ptl1_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:
	@echo "VALID"
