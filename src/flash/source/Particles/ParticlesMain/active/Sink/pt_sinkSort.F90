!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkSort
!!
!! NAME
!!
!!  pt_sinkSort
!!
!! SYNOPSIS
!!
!!  call NewQsort_IN(id_sorted, QSindex)
!!
!! DESCRIPTION
!!
!!  Quick sort routine to sort sink particles by tag. Returns the indexes
!!  of the original array, sorted by tag (QSindex). This is used before
!!  reducing particle properties across all processors, such that each
!!  processor reduces an ordered list of particles (instead of the un-ordered
!!  particles_global list).
!!
!! ARGUMENTS
!!
!!   id_sorted - the sorted particle indexes
!!
!!   QSindex - the indexes of the sorted particles
!!
!! NOTES
!!
!!  COPIED AND MODIFIED 2010 for sink particles in FLASH!
!!
!!  A MODULE for non numerical routines (sorting and locate)
!!
!!  Copyright (C) 2005  Alberto Ramos <alberto@martin.ft.uam.es>
!!
!!  This program is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 2 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program; if not, write to the Free Software
!!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
!!
!!***

MODULE pt_sinkSort
! *
! ***************************************************
! *
! * Quicksort routines for sorting 
! *
! ***************************************************


  IMPLICIT NONE

  Private Insrt_IN, Swap_IN, Partition_IN

CONTAINS


! ***********************************
! *
  Subroutine Insrt_IN(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. Integer 
! * version.
! ***********************************
    implicit none
    
    Integer(8), Intent (inout) :: X(:)
    Integer(8), Intent (out), Optional :: Ipt(:)
    
    Integer(8) :: Rtmp
    Integer(8) :: I, J


    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
       
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL Swap_IN(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine Insrt_IN




! ***********************************
! *
  Subroutine Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
    implicit none
    
    Integer(8), Intent (inout) :: X(:)
    Integer(8), Intent (in) :: I, J

    Integer(8) :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap_IN


! ***********************************
! *
  Subroutine NewQsort_IN(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt.
! ***********************************
    implicit none
    
    Type Limits
       Integer(8) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(8), Parameter :: Isw = 10

    Integer(8), Intent (inout) :: X(:)
    Integer(8), Intent (out), Optional :: Ipt(:)
    
    Integer(8) :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Type (Limits), Allocatable :: Stack(:)
    
    
    Allocate(Stack(2*Size(X)))



    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition_IN(X, Ileft, Iright, Ipvn, Ipt)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
       
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrt_IN(X(Ileft:Iright))
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition_IN(X, Ileft, Iright, Ipvn)
             
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    End If

    Deallocate(Stack)

    Return
    
  CONTAINS

    ! ***********************************
    Integer(8) Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for pt_sinkSort.
    ! ***********************************
      implicit none
      
      Integer(8), Intent (in) :: XX(:)
      Integer(8), Intent (in) :: IIleft, IIright

!       Integer(8), Parameter :: DP = Kind(1.D0)
      Integer(8) :: XXcp(3)    ! used to be    Real (kind=DP)
      Integer(8) :: IIpt(3), IImd
      
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IIright)
      XXcp(3) = XX(IImd)
      
      CALL Insrt_IN(XXcp, IIpt)
      
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return
    End Function ChoosePiv

    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************
      implicit none
      
      Integer(8), Intent (inout) :: XX(:)
      Integer(8), Intent (inout) :: IIpt(:)
      Integer(8), Intent (in) :: IIl, IIr
      
      Integer(8) :: RRtmp
      Integer(8) :: II, JJ
      

      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap_IN(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
      
      Return
    End Subroutine InsrtLC


  End Subroutine NewQsort_IN


! ***********************************
! *
  Integer(8) Function Partition_IN(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by pt_sinkSort.
! ***********************************
    implicit none
    
    Integer(8), Intent (inout) :: X(:)
    Integer(8), Intent (in) :: Ileft, Iright, Ipv
    Integer(8), Intent (inout), Optional :: Ipt(:)
    
    Integer(8) :: Rpv
    Integer(8) :: I

    Rpv = X(Ipv)
    CALL Swap_IN(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap_IN(X, I, Ipvfn)
             CALL Swap_IN(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap_IN(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL Swap_IN(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)

    Return
  End Function Partition_IN


End MODULE pt_sinkSort

