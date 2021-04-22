#!/usr/bin/perl -w
use SVG;
use strict;
use POSIX;
$|=1;

my $DRAW_TICS = 1;
my $xBin = 100;
my $yBin = 100;
my $textY = 6.00;
my $delta = 0.15;
#my $FDR = 0.05;
my $FDR = 1e-4;
my $motifWidth = 8;
my $usage = "Usage:\n$0 <x y z data file> [max delta]\n";

my $mapFile = $ARGV[0] or die $usage;
#my $mode = defined($ARGV[1]) ? $ARGV[1] : die $usage;
my $mode = 0;
my $maxDelta = defined($ARGV[1]) ? $ARGV[1] : 0;
my $pValueMapFile = $ARGV[2] ? $ARGV[2] : "";
print "building heat maps for $mapFile\n";

my($maxZ,$data) = readMapFile($mapFile,$maxDelta);
my($width,$height,$xOffset) = getDataRange($data);
buildScaleBar();

my $fontHash = {'font' => 'Arial',
		'font-size' => 70};
my $image = SVG->new(width => ($width+1)*$xBin, height => ($height+2+$textY)*$yBin);
my $COUNT = 0;

foreach my $datum (@{$data}) {
    my($x,$y,$z) = @{$datum};
    my $xPos = ($x + $xOffset)*$xBin;
    my $yPos = $y*$yBin;
    my($r,$g,$b) = getColor($z,$maxZ); 
    $image->rectangle('x' => $xPos, 'y' => $yPos, 'width' => $xBin, 'height' => $yBin, style => {'fill' => "rgb($r,$g,$b)" , 'stroke' => "rgb(255,255,255)", 'stroke-width' => 0 });
    if($DRAW_TICS) {
	if($y == 4) {
	    $image->line('x1' => $xPos + $xBin/2, 'y1' => $yPos + $yBin, 'x2' => $xPos + $xBin/2, 'y2' => $yPos + 1.25*$yBin, 'stroke' => "rgb(0,0,0)", 'stroke-width' => 5);
	}
    }
}

if($DRAW_TICS) {
    for(my $i=-$xOffset;$i<=$xOffset+1;$i++) {
	my $xPos = ($xOffset+$i+0.25)*$xBin;
	$xPos = ($xOffset+$i)*$xBin if($i <= -1);
	$xPos = ($xOffset+$i-0.25)*$xBin if($i <= -10);
	$xPos = ($xOffset+$i)*$xBin if($i >= 10);
	#my $pos = $i;
	my $pos = ($i <= 0) ? $i - 1 : $i;
	#print $pos . "\t". $xOffset ."\n";
	#if($pos % 3 == 0) {
	$image->text('style' => $fontHash, 'x' => $xPos, 'y' => $textY*$yBin,  -cdata => $pos);
	#}
    }
}

for(my $i=-$xOffset;$i<=$xOffset;$i++) {
    if($mode =~ /2/) {
	if($i % 3 == 0 and $i <= 0) {
	    my $xPos1 = ($xOffset+$i+$delta)*$xBin;    
	    my $xPos2 = ($xOffset+$i+3-$delta)*$xBin;
	    my $yPos = 4.2*$yBin;
	    $image->line('x1' => $xPos1, 'y1' => $yPos, 'x2' => $xPos2, 'y2' => $yPos, 'stroke' => "rgb(0,0,0)", 'stroke-width' => 10);
	}    	
    } 
    if($mode =~ /1/) {
	if($i % 3 == 0 and $i >= 0) {
	    my $xPos1 = ($xOffset+$i+$delta)*$xBin;    
	    my $xPos2 = ($xOffset+$i+3-$delta)*$xBin;
	    my $yPos = 4.2*$yBin;
	    $image->line('x1' => $xPos1, 'y1' => $yPos, 'x2' => $xPos2, 'y2' => $yPos, 'stroke' => "rgb(0,0,0)", 'stroke-width' => 10);
	}    
    }
    if($mode =~ /3/) {
	if($i == 0) {
	    my $xPos1 = ($xOffset+$i+$delta)*$xBin;    
	    my $xPos2 = ($xOffset+$i+$motifWidth-$delta)*$xBin;
	    my $yPos = 4.2*$yBin;
	    $image->line('x1' => $xPos1, 'y1' => $yPos, 'x2' => $xPos2, 'y2' => $yPos, 'stroke' => "rgb(0,0,0)", 'stroke-width' => 10);	    
	}
    }
}

if($pValueMapFile) {
    my $xShift = 0.25;
    my $pValueMap = readPValueMapFile($pValueMapFile);
    for(my $i=-$xOffset;$i<=$xOffset;$i++) {
	for(my $j=0;$j<5;$j++) {
	    my $xPos = ($i + $xOffset + $xShift)*$xBin;
	    my $yPos = ($j + 1)*$yBin;
	    if($pValueMap->{$i}{$j}) {
		$image->text('style' => $fontHash, 'x' => $xPos, 'y' => $yPos, -cdata => "*");
	    }
	}
    }
}


my($out) = $mapFile =~ /(.*)\.txt/;
#my $outputFile = $pValueMapFile ? $out . ".qVal.svg" : $out . ".svg";
my $outputFile = $out . ".svg";
open(OUTPUT,">$outputFile");
print OUTPUT $image->xmlify();
close(OUTPUT);

sub buildScaleBar {
    my $scaleWidth = 50;
    my $scaleBuffer = 50;
    my $scaleHeight = 300;
    my $scaleYBuffer = 15;
    my $scale = SVG->new(width => $scaleWidth+2*$scaleBuffer, height => $scaleHeight+4*$scaleYBuffer);
    my $COUNT = 0;
    my $deltaZ = 2*$maxZ/10;
    my $i=0;
    for(my $Z = 0; $Z <= 2*$maxZ+0.001; $Z += $deltaZ) {	
	my $xPos = $scaleBuffer;
	my $yPos = 2*$scaleYBuffer + $scaleHeight - ($scaleYBuffer + ($Z/(2*$maxZ))*$scaleHeight);
	my $Zval = $Z - $maxZ;
	
	my($r,$g,$b) = getColor($Zval,$maxZ);
	$scale->rectangle('x' => $xPos, 'y' => $yPos, 'width' => $scaleWidth, 'height' => ($deltaZ/(2*$maxZ))*$scaleHeight, style => {'fill' => "rgb($r,$g,$b)" , 'stroke' => "rgb(255,255,255)", 'stroke-width' => 0 });
	my $textYPos = 2*$scaleYBuffer + $scaleHeight - ($scaleYBuffer + ($Z/(2*$maxZ))*$scaleHeight) + 0.7*($scaleHeight/11);
	$scale->text('x' => $scaleWidth+$scaleBuffer, 'y' => $textYPos, -cdata => sprintf("%.2f",$Zval));
	$i++;
    }
    my($out) = $mapFile =~ /(.*)\.txt/;
    #my $outputFile = $pValueMapFile ? $out . ".scalebar.qVal.svg" : $out . ".scalebar.svg";
    my $outputFile = $out . ".scalebar.svg";
    open(OUTPUT,">$outputFile");
    print OUTPUT $scale->xmlify();
    close(OUTPUT);
}

sub readPValueMapFile {
    my($pValueMapFile) = @_;
    my %sig;
    open(PVALUE,$pValueMapFile) or die "Could not open $pValueMapFile";
    while(<PVALUE>) {
	chomp;
	my($pos,$base,$pValue,$qValue) = split(/\t/);
	if($qValue <= $FDR) {
	    $sig{$pos}{$base} = 1;
	}
    }
    return \%sig;
}

sub getColor {
    my($z,$maxZ) = @_;
    my $r = 255;
    my $g = 255;
    my $b = 255;
    if($z > 0) {
	# positive scores get red
	$g = int 255*(($maxZ-abs($z))/$maxZ);
	$b = int 255*(($maxZ-abs($z))/$maxZ);
    } else {
	# negative scores get blue
	$r = int 255*(($maxZ-abs($z))/$maxZ);
	$g = int 255*(($maxZ-abs($z))/$maxZ);
    }
    return ($r,$g,$b);
}

sub getDataRange {
    my($data) = @_;
    my $minX = 10e10;
    my $maxX = -10e10;
    my $minY = 10e10;
    my $maxY = -10e10;
    my $minZ = 10e10;
    my $maxZ = -10e10;
    foreach my $datum (@{$data}) {
	my($x,$y,$z) = @{$datum};
	$maxX = $x if($x > $maxX);
	$minX = $x if($x < $minX);
	$maxY = $y if($y > $maxY);
	$minY = $y if($y < $minY);
	$maxZ = $z if($z > $maxZ);
	$minZ = $z if($z < $minZ);
    }
    return $maxX-$minX, $maxY-$minY, abs($minX)
}

sub max {
    my($a,$b) = @_;
    if($a > $b) {
	return $a;
    }
    return $b;
}

sub readMapFile {
    my($mapFile,$maxDelta) = @_;
    my $maxZ;
    my @data;
    open(FILE,$mapFile) or die "could not open $mapFile\n";
    while(<FILE>) {
	chomp;
	if(/^#/) {
	    if(/^#max/) {
		($maxZ) = $_ =~ /max=(.*)/;
	    }
	} else {
	    my($x,$y,$z) = split(/\t/);
	    push(@data,[$x,$y,$z]);
	}
    }
    if ($maxDelta) {
	$maxZ = $maxDelta;
    }
    return($maxZ,\@data);
}
