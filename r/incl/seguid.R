## SEGUID v1 on linear single-stranded DNA
seguid("GATTACA")
#> seguid=tp2jzeCM2e3W4yxtrrx09CMKa/8

## SEGUID v2 on linear single-stranded DNA
lsseguid("GATTACA")
#> lsseguid=tp2jzeCM2e3W4yxtrrx09CMKa_8

## SEGUID v2 on cicular single-stranded DNA
## GATTACA = ATTACAG = ... = AGATTAC
csseguid("GATTACA")
#> csseguid=mtrvbtuwr6_MoBxvtm4BEpv-jKQ

## SEGUID v2 on blunt, linear double-stranded DNA
##   GATTACA
##   CTAATGT
ldseguid("GATTACA", "TGTAATC", overhang = 0)
#> ldseguid=AcRsEcNFrui5wCxI7xxo6wnDYPY

## SEGUID v2 on staggered, linear double-stranded DNA
##   -ATTACA
##   CTAAT--
ldseguid("-ATTACA", "--TAATC")
#> ldseguid=98Klwxd3ZQPGHqnH3BheIuZVHQQ

## SEGUID v2 on circular double-stranded DNA
## GATTACA = ATTACAG = ... = AGATTAC
## CTAATGT = TAATGTC = ... = TCTAATG
cdseguid("GATTACA", "TGTAATC")
#> cdseguid=zCuq031K3_-40pArbl-Y4N9RLnA

## SEGUID v2 on linear single-stranded expanded
## epigenetic sequence (Viner et al., 2024)
viner_DNA <- "{DNA},m1,1m,h2,2h,f3,3f,c4,4c"
lsseguid("AmT2C", alphabet = viner_DNA)
#> lsseguid=MW4Rh3lGY2mhwteaSKh1-Kn2fGA

## SEGUID v2 on linear double-stranded expanded
## epigenetic sequence (Viner et al., 2024)
ldseguid("AmT2C", "GhA1T", overhang = 0, alphabet = viner_DNA)
#> ldseguid=rsPDjP4SWr3-ploCeXTdTA80u0Y
