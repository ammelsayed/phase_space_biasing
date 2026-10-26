C ************************************************************
C Source for the library implementing a bias function that
C biases the cross section itself instead of a particular
C kinematic distribution: it flattens the event generation
C rate over the whole phase space, so that regions with a
C small differential cross section are populated with as many
C events as regions with a large one (relative oversampling
C of the low-cross-section regions).
C
C It uses the original weight of the event, which is
C proportional to the differential cross section at that
C phase-space point:
C    bias_wgt = xsec_flat_bias_reference_weight
C               / |original_weight|^xsec_flat_bias_exponent
C so that the total weight
C    original_weight * bias_weight
C is (approximately) constant for exponent=1. The generator
C then samples the underlying integration measure uniformly
C and the unweighting efficiency approaches unity.
C
C Values 0 < xsec_flat_bias_exponent < 1 only partially
C flatten the sampling; exponent=0 recovers the unbiased
C generation.
C
C This module is used with impact_xsec=.False. (set below) so
C that the physical cross section and distributions are NOT
C modified: the unweighted events are simply distributed
C differently over the phase space, and the bias weight of
C each event is stored in the <rwgt> block of the LHE file so
C that any analysis can undo (or change) the bias a
C posteriori. The total cross section reported by MG5aMC
C remains the unbiased one.
C
C /!\ Experimental module: the flattening is performed in the
C integration measure of the generator (which already includes
C its own importance sampling), NOT uniformly in a physical
C observable. Since the original weight can span many orders
C of magnitude, the bias weight can become very large; if the
C unweighting efficiency degrades, reduce the exponent.
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the
C relevant parameters for this bias module.
C
C  parameters = {'xsec_flat_bias_exponent': 1.0,
C               'xsec_flat_bias_reference_weight': 1.0}
C

      subroutine bias_wgt(p, original_weight, bias_weight)
          implicit none
C
C Parameters
C
          include '../../maxparticles.inc'
          include '../../nexternal.inc'

C
C Arguments
C
          double precision p(0:3,nexternal)
          double precision original_weight, bias_weight
C
C local variables
C
c
c local variables defined in the run_card
c
          double precision xsec_flat_bias_exponent
          double precision xsec_flat_bias_reference_weight
C
C Global variables
C
C
C Mandatory common block to be defined in bias modules
C
          double precision stored_bias_weight
          data stored_bias_weight/1.0d0/
          logical impact_xsec, requires_full_event_info
C         We do NOT want to modify the physical cross section: the
C         bias is only used to redistribute the generated events
C         over the phase space, and it is recorded event by event
C         in the <rwgt> block of the LHE file.
          data impact_xsec/.False./
C         Of course this module does not require the full event
C         information (color, resonances, helicities, etc..)
          data requires_full_event_info/.False./
          common/bias/stored_bias_weight,impact_xsec,
     &                requires_full_event_info

C
C    Setup the value of the parameters from the run_card
C
      include '../bias.inc'

C --------------------
C BEGIN IMPLEMENTATION
C --------------------

C         Protect against vanishing original weights
          if (abs(original_weight).lt.1.0d-99) then
            bias_weight = 1.0d0
            return
          endif

          bias_weight = xsec_flat_bias_reference_weight
     &         / abs(original_weight)**xsec_flat_bias_exponent

          return

      end subroutine bias_wgt
