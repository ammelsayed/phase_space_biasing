include ../../make_opts

all: ht_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module ht_bias
#
ht_bias.o: ht_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o ht_bias.o ht_bias.f

ht_bias: ht_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:
	@echo "VALID"
