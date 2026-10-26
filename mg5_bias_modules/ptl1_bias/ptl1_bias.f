C ************************************************************
C Source for the library implementing a bias function that
C populates the large pt tail of the leading charged lepton.
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) ptl1_bias_target_ptl1 : target ptl1 value
C    > (double precision) ptl1_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [ptl1(evt)/ptl1_bias_target_ptl1]^enhancement_power
C
C where ptl1(evt) is the transverse momentum of the hardest
C charged lepton (e+-, mu+-, ta+-) of the event. Events
C without any charged lepton are left unbiased.
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the
C relevant parameters for this bias module.
C
C  parameters = {'ptl1_bias_target_ptl1': 300.0,
C               'ptl1_bias_enhancement_power': 4.0}
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
c
c local variables defined in the run_card
c
          double precision ptl1_bias_target_ptl1
          double precision ptl1_bias_enhancement_power
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

          do i=1,nexternal
            ptl1(i)=-1.0d0
            if (is_a_l(i)) then
              ptl1(i)=sqrt(p(1,i)**2+p(2,i)**2)
            endif
          enddo

          max_ptl1=-1.0d0
          do i=1,nexternal
            max_ptl1 = max(max_ptl1,ptl1(i))
          enddo
          if (max_ptl1.lt.0.0d0) then
            bias_weight = 1.0d0
            return
          endif

          bias_weight = (max_ptl1/ptl1_bias_target_ptl1)
     &                                   **ptl1_bias_enhancement_power

          return

      end subroutine bias_wgt
