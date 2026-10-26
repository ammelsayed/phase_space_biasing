C ************************************************************
C Source for the library implementing a bias function that
C populates simultaneously the large tails of the transverse
C momentum of the leading charged lepton (ptl1), of the
C missing transverse energy (met) and of the scalar sum of the
C transverse momenta of all visible particles (ht). It is the
C product of the three individual bias modules ptl1_bias,
C met_bias and ht_bias.
C
C The options of this subroutine, that can be set in
C the run card are:
C    > (double precision) ptl1_bias_target_ptl1 : target ptl1 value
C    > (double precision) ptl1_bias_enhancement_power : ptl1 exponent
C    > (double precision) met_bias_target_met : target met value
C    > (double precision) met_bias_enhancement_power : met exponent
C    > (double precision) ht_bias_target_ht : target ht value
C    > (double precision) ht_bias_enhancement_power : ht exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [ptl1(evt)/ptl1_bias_target_ptl1]^ptl1_power
C             * [met(evt)/met_bias_target_met]^met_power
C             * [ht(evt)/ht_bias_target_ht]^ht_power
C
C where ptl1(evt) is the pt of the hardest charged lepton,
C met(evt) is the missing transverse energy (magnitude of the
C vector sum of the transverse momenta of the invisible
C particles) and ht(evt) is the scalar sum of the transverse
C momenta of all the visible final-state particles. A factor
C is set to unity when its observable is not defined for the
C event (no charged lepton / no invisible particle / no
C visible transverse momentum).
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the
C relevant parameters for this bias module.
C
C  parameters = {'ptl1_bias_target_ptl1': 300.0,
C               'ptl1_bias_enhancement_power': 4.0,
C               'met_bias_target_met': 300.0,
C               'met_bias_enhancement_power': 4.0,
C               'ht_bias_target_ht': 1000.0,
C               'ht_bias_enhancement_power': 4.0}
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
          double precision ptl1(nexternal)
          double precision max_ptl1
          double precision metx, mety, met
          double precision ht
          double precision wgt_ptl1, wgt_met, wgt_ht
c
c local variables defined in the run_card
c
          double precision ptl1_bias_target_ptl1
          double precision ptl1_bias_enhancement_power
          double precision met_bias_target_met
          double precision met_bias_enhancement_power
          double precision ht_bias_target_ht
          double precision ht_bias_enhancement_power
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

C         Bias on the pt of the leading charged lepton
          max_ptl1=-1.0d0
          do i=1,nexternal
            ptl1(i)=-1.0d0
            if (is_a_l(i)) then
              ptl1(i)=sqrt(p(1,i)**2+p(2,i)**2)
              max_ptl1 = max(max_ptl1,ptl1(i))
            endif
          enddo
          if (max_ptl1.lt.0.0d0) then
            wgt_ptl1 = 1.0d0
          else
            wgt_ptl1 = (max_ptl1/ptl1_bias_target_ptl1)
     &                                   **ptl1_bias_enhancement_power
          endif

C         Bias on the missing transverse energy
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
            wgt_met = 1.0d0
          else
            wgt_met = (met/met_bias_target_met)
     &                                   **met_bias_enhancement_power
          endif

C         Bias on the scalar sum of the visible transverse momenta
          ht = 0.0d0
          do i=1,nexternal
            if (.not.is_a_nu(i)) then
              ht = ht + sqrt(p(1,i)**2+p(2,i)**2)
            endif
          enddo
          if (ht.le.0.0d0) then
            wgt_ht = 1.0d0
          else
            wgt_ht = (ht/ht_bias_target_ht)
     &                                    **ht_bias_enhancement_power
          endif

          bias_weight = wgt_ptl1*wgt_met*wgt_ht

          return

      end subroutine bias_wgt
