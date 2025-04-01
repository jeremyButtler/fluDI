#!/usr/bin/env Rscript

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec01:
#   - variable declarations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inputStr = commandArgs(); # Holds the user input
fileStr = "2024-mooflu-SRR29165088-map-filt-tbCon-blast-fluDI-diRna-diCoords.tsv";

# flu segment lengths
pb2LenSI = 2341;
pb1LenSI = 2341;
paLenSI = 2233;
haLenSI = 1778;
npLenSI = 1565;
naLenSI = 1413;
mLenSI = 1027;
nsLenSI = 890;

pb2Bl = 0;
pb1Bl = 0;
paBl = 0;
haBl = 0;
npBl = 0;
naBl = 0;
mBl = 0;
nsBl = 0;

pb2DelArySI = rep(0, pb2LenSI);
pb1DelArySI = rep(0, pb1LenSI);
paDelArySI = rep(0, paLenSI);
haDelArySI = rep(0, haLenSI);
npDelArySI = rep(0, npLenSI);
naDelArySI = rep(0, naLenSI);
mDelArySI = rep(0, mLenSI);
nsDelArySI = rep(0, nsLenSI);

pb2CoverArySI = rep(0, pb2LenSI);
pb1CoverArySI = rep(0, pb1LenSI);
paCoverArySI = rep(0, paLenSI);
haCoverArySI = rep(0, haLenSI);
npCoverArySI = rep(0, npLenSI);
naCoverArySI = rep(0, naLenSI);
mCoverArySI = rep(0, mLenSI);
nsCoverArySI = rep(0, nsLenSI);

zeroArySI = rep(0, pb1LenSI);
depthSI = 0;
skipDelSI = 20; # bases skipped at ends
diArySI = NULL; # for deletion depth

outNameStr = "out";
widthSI = 1400;
heightSI = 500;

if(length(inputStr) < 6){
  print("nothing input");
  print("Rscript diGraph.r DI-coords.tsv prefix");

} else if(length(inputStr) < 7) {
  print("not diCoords tsv file input");
  print("Rscript diGraph.r prefix DI-coords.tsv");

} else { # Else: have input

outNameStr = inputStr[6];
fileStr = inputStr[7];

if(! file.exists(fileStr) ) {
   print( paste("could not open", fileStr) );
} else { # Else: could read in file

diDF =
  read.csv(
     fileStr,
     sep = "\t",
     header = TRUE
   );

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec02:
#   - find depths (assuming complete mapping)
#   o main sec02 sub01:
#     - get deletion depth
#   o main sec02 sub02:
#     - get coverage depth
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*********************************************************
# Main Sec02 Sub01:
#   - get deletion depth
#*********************************************************

for( siRow in 1:length(diDF$read) )
{ # Loop: find number reads supporting a deletion

   diArySI = diDF$start_di[siRow]:diDF$end_di[siRow];

   if( is.na(diDF$ref[siRow]) ) {
      naBl = 1;
      naDelArySI[diArySI] = naDelArySI[diArySI] + 1;
      # in this case "NA" is converted to NA

   } else if(diDF$ref[siRow] == "PB2") {
      pb1Bl = 1;
      pb1DelArySI[diArySI] = pb1DelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "PB1") {
      pb2Bl = 1;
      pb2DelArySI[diArySI] = pb2DelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "PA") {
      paBl = 1;
      paDelArySI[diArySI] = paDelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "HA") {
      haBl = 1;
      haDelArySI[diArySI] = haDelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "NP") {
      npBl = 1;
      npDelArySI[diArySI] = npDelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "NA") {
      naBl = 1;
      naDelArySI[diArySI] = naDelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "M") {
      mBl = 1;
      mDelArySI[diArySI] = mDelArySI[diArySI] + 1;

   } else if(diDF$ref[siRow] == "NS") {
      nsBl = 1;
      nsDelArySI[diArySI] = nsDelArySI[diArySI] + 1;

   } # find segment

} # Loop: find number reads supporting a deletion

#*********************************************************
# Main Sec02 Sub02:
#   - get coverage depth
#*********************************************************

coverDF = diDF[! duplicated(diDF$read), ];# get unique ids

coverDF;

for( siRow in 1:length(coverDF$read) )
{ # Loop: get read depth

   covArySI =
      coverDF$ref_start[siRow]:coverDF$ref_end[siRow];

   if( is.na(coverDF$ref[siRow]) ) {
      naCoverArySI[covArySI] = naCoverArySI[covArySI] + 1;
      # in this case "NA" is converted to NA

   } else if(coverDF$ref[siRow] == "PB2") {
      pb1CoverArySI[covArySI]= pb1CoverArySI[covArySI] +1;

   } else if(coverDF$ref[siRow] == "PB1") {
      pb2CoverArySI[covArySI]= pb2CoverArySI[covArySI] +1;

   } else if(coverDF$ref[siRow] == "PA") {
      paCoverArySI[covArySI] = paCoverArySI[covArySI] + 1;

   } else if(coverDF$ref[siRow] == "HA") {
      haCoverArySI[covArySI] = haCoverArySI[covArySI] + 1;

   } else if(coverDF$ref[siRow] == "NP") {
      npCoverArySI[covArySI] = npCoverArySI[covArySI] + 1;

   } else if(coverDF$ref[siRow] == "NA") {
      naCoverArySI[covArySI] = naCoverArySI[covArySI] + 1;

   } else if(coverDF$ref[siRow] == "M") {
      mCoverArySI[covArySI] = mCoverArySI[covArySI] + 1;

   } else if(coverDF$ref[siRow] == "NS") {
      nsCoverArySI[covArySI] = nsCoverArySI[covArySI] + 1;

   } # find segment

} # Loop: get read depth

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec03:
#   - PB1 graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(pb1Bl != 0)
{ # If: have pb1 segment
   png(
      filename = paste(outNameStr, "PB1.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(pb1CoverArySI),
      y = pb1CoverArySI,
      type = "l",                         # line graph
      xlab = "PB1 position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:pb1LenSI,                     # positions
      y = pb1DelArySI                     # read depths
   ); # plot coverage

   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(pb1DelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(pb1LenSI - skipDelSI, pb1LenSI - skipDelSI),
      y = c( 0, max(pb1DelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have pb1 segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec04:
#   - PB2 graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(pb2Bl != 0)
{ # If: have pb2 segment
   png(
      filename = paste(outNameStr, "PB2.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(pb2CoverArySI),
      y = pb2CoverArySI,
      type = "l",                         # line graph
      xlab = "PB2 position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:pb2LenSI,                     # positions
      y = pb2DelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(pb2DelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(pb2LenSI - skipDelSI, pb2LenSI - skipDelSI),
      y = c( 0, max(pb2DelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have pb2 segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec05:
#   - PA graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(paBl != 0)
{ # If: have pa segment
   png(
      filename = paste(outNameStr, "PA.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(paCoverArySI),
      y = paCoverArySI,
      type = "l",                         # line graph
      xlab = "PA position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:paLenSI,                     # positions
      y = paDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(paDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(paLenSI - skipDelSI, paLenSI - skipDelSI),
      y = c( 0, max(paDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have pa segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec06:
#   - HA graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(haBl != 0)
{ # If: have ha segment
   png(
      filename = paste(outNameStr, "HA.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(haCoverArySI),
      y = haCoverArySI,
      type = "l",                         # line graph
      xlab = "HA position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:haLenSI,                     # positions
      y = haDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(haDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(haLenSI - skipDelSI, haLenSI - skipDelSI),
      y = c( 0, max(haDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have ha segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec07:
#   - NP graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(npBl != 0)
{ # If: have np segment
   png(
      filename = paste(outNameStr, "NP.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(npCoverArySI),
      y = npCoverArySI,
      type = "l",                         # line graph
      xlab = "NP position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:npLenSI,                     # positions
      y = npDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(npDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(npLenSI - skipDelSI, npLenSI - skipDelSI),
      y = c( 0, max(npDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have np segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec09:
#   - NA graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(naBl != 0)
{ # If: have na segment
   png(
      filename = paste(outNameStr, "NA.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(naCoverArySI),
      y = naCoverArySI,
      type = "l",                         # line graph
      xlab = "NA position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:naLenSI,                     # positions
      y = naDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(naDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(naLenSI - skipDelSI, naLenSI - skipDelSI),
      y = c( 0, max(naDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have na segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec10:
#   - M graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(mBl != 0)
{ # If: have m segment
   png(
      filename = paste(outNameStr, "M.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(mCoverArySI),
      y = mCoverArySI,
      type = "l",                         # line graph
      xlab = "M position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:mLenSI,                     # positions
      y = mDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(mDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(mLenSI - skipDelSI, mLenSI - skipDelSI),
      y = c( 0, max(mDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have m segment

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Sec11:
#   - NS graph
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(nsBl != 0)
{ # If: have ns segment
   png(
      filename = paste(outNameStr, "NS.png", sep = "-"),
      widthSI,
      heightSI
   ); # file to save graph as
   
   plot(
      x = 1:length(nsCoverArySI),
      y = nsCoverArySI,
      type = "l",                         # line graph
      xlab = "NS position",              # x-axis label
      ylab = "Number reads with deletion",# y-axis label
      col = "DARKGREY"
   ); # maximum read depth plot
   
   lines(
      x = 1:nsLenSI,                     # positions
      y = nsDelArySI                     # read depths
   ); # plot coverage
   
   # regions ignored large deletions in
   lines(
      x = c(skipDelSI, skipDelSI),
      y = c( 0, max(nsDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
   
   lines(
      x = c(nsLenSI - skipDelSI, nsLenSI - skipDelSI),
      y = c( 0, max(nsDelArySI) ),
      lty = 2,
      lwd = 1.5
   );
} # If: have ns segment

} # Else: could read in file
} # Else: have input
