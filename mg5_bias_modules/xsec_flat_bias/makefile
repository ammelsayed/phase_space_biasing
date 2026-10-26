include ../../make_opts

all: xsec_flat_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module xsec_flat_bias
#
xsec_flat_bias.o: xsec_flat_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o xsec_flat_bias.o xsec_flat_bias.f

xsec_flat_bias: xsec_flat_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:
	@echo "VALID"
