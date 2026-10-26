C ************************************************************
C Source for the library implementing a bias function that
C populates the large tail of the scalar sum of the transverse
C momenta of all the visible particles of the event (HT).
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) ht_bias_target_ht : target ht value
C    > (double precision) ht_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [ht(evt)/ht_bias_target_ht]^enhancement_power
C
C where ht(evt) is the scalar sum of the transverse momenta of
C all the visible final-state particles (jets, b-jets, charged
C leptons and photons), i.e. all particles except the
C invisible ones (neutrinos). Events with no visible
C transverse momentum are left unbiased.
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the
C relevant parameters for this bias module.
C
C  parameters = {'ht_bias_target_ht': 1000.0,
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
          double precision ht
c
c local variables defined in the run_card
c
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

          ht = 0.0d0
          do i=1,nexternal
            if (.not.is_a_nu(i)) then
              ht = ht + sqrt(p(1,i)**2+p(2,i)**2)
            endif
          enddo

          if (ht.le.0.0d0) then
            bias_weight = 1.0d0
            return
          endif

          bias_weight = (ht/ht_bias_target_ht)
     &                                    **ht_bias_enhancement_power

          return

      end subroutine bias_wgt
