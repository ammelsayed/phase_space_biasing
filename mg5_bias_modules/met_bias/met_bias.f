C ************************************************************
C Source for the library implementing a bias function that
C populates the large tail of the missing transverse energy.
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) met_bias_target_met : target met value
C    > (double precision) met_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [met(evt)/met_bias_target_met]^enhancement_power
C
C where met(evt) is the missing transverse energy of the
C event, computed as the magnitude of the vector sum of the
C transverse momenta of the invisible particles (neutrinos).
C Events without any invisible particle are left unbiased.
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the
C relevant parameters for this bias module.
C
C  parameters = {'met_bias_target_met': 300.0,
C               'met_bias_enhancement_power': 4.0}
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
          integer i
          double precision metx, mety, met
c
c local variables defined in the run_card
c
          double precision met_bias_target_met
          double precision met_bias_enhancement_power
C
C Global variables
C
C
C Mandatory common block to be defined in bias modules
C
          double precision stored_bias_weight
          data stored_bias_weight/1.0d0/
          logical impact_xsec, requires_full_event_info
C         We only want to bias distributions, but not impact the xsec.
          data impact_xsec/.False./
C         Of course this module does not require the full event
C         information (color, resonances, helicities, etc..)
          data requires_full_event_info/.False./
          common/bias/stored_bias_weight,impact_xsec,
     &                requires_full_event_info
C
C Accessingt the details of the event
C
          logical is_a_j(nexternal),is_a_l(nexternal),
     &            is_a_b(nexternal),is_a_a(nexternal),
     &            is_a_onium(nexternal),is_a_nu(nexternal),
     &            is_heavy(nexternal),do_cuts(nexternal)
          common/to_specisa/is_a_j,is_a_a,is_a_l,is_a_b,is_a_nu,
     &                      is_heavy,is_a_onium,do_cuts

C
C    Setup the value of the parameters from the run_card
C
      include '../bias.inc'

C --------------------
C BEGIN IMPLEMENTATION
C --------------------

          metx = 0.0d0
          mety = 0.0d0
          do i=1,nexternal
            if (is_a_nu(i)) then
              metx = metx - p(1,i)
              mety = mety - p(2,i)
            endif
          enddo
          met = sqrt(metx**2+mety**2)

          if (met.le.0.0d0) then
            bias_weight = 1.0d0
            return
          endif

          bias_weight = (met/met_bias_target_met)
     &                                   **met_bias_enhancement_power

          return

      end subroutine bias_wgt
