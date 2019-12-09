#!/usr/bin/perl -w

############################################################
#
#  To Run
#  for i in {1..500};do echo $i; perl MIF2CLconversion-V2.pl Cathepsin N:-6,O:-2.7,OH:-5.5 CathD-model-$i\_ion.kont;done
#
#
###########################################################



use warnings;
use strict; 

my$target_name = $ARGV[0];												## get the protein target name
chomp($target_name);

my$Probe_Type = $ARGV[1];
chomp($Probe_Type);
my@Probe_Type_array = split(",",$Probe_Type);
my$Probe_type_length = scalar(@Probe_Type_array);

my$Filename = $ARGV[2];
chomp($Filename);

open(A,$Filename);
my@MIF = <A>;

for(my$i=0;$i<scalar@Probe_Type_array;$i++)
{
	my$r = $i+1;
	my($x,$y,$z,$e) = probe_XYZE(\@MIF,$Probe_type_length,$r);
	#XYZE2mol2($x,$y,$z,$e,$Probe_Type_array[$i],$Filename);
	CliquePharm_Input($x,$y,$z,$e,$Probe_Type_array[$i],$target_name,$Filename);
}



sub probe_XYZE
{
	my($MIF_data,$NP,$Probe_order) = @_;
	my@X = ();
	my@Y = ();
	my@Z = ();
	my@E = ();
	
	
	my$Grid_PPP = (scalar(@{$MIF_data}))/$NP;
	my$start_point = $Grid_PPP*($Probe_order-1);
	my$End_point   = $Grid_PPP*$Probe_order;
	
	my$MIF_Plane1 = int(($End_point - $start_point)/2);
	my$MIF_Plane2 = $start_point + $MIF_Plane1;
	
	for(my$i = $start_point;$i<$MIF_Plane2;$i++)
	{
		my@tx = split(" ",$$MIF_data[$i]);
		push(@X,$tx[1]);
		push(@Y,$tx[2]);
		push(@Z,$tx[3]);
	}
	for(my$j = $MIF_Plane2+1;$j < $End_point;$j++)
	{
		my@te = split(" ",$$MIF_data[$j]);
		push(@E,$te[0]);
	}
	return(\@X,\@Y,\@Z,\@E);
}

	
sub XYZE2mol2
{
	my($X,$Y,$Z,$E,$Probe_Type,$name) = @_;
	my@S = split(":",$Probe_Type);
	my$Number_Points = scalar(@$X);
	my@Tmp = split("_",$name);
	print "$Tmp[0]\t$S[0]\n";  
	my$file = "$Tmp[0]_$S[0].mol2";
	open(OUT,">../CliqueP_CathD/$file");                            ### Change Output Directory Name.
	print OUT "@<TRIPOS>MOLECULE\n";
	print OUT "$Number_Points 0 1 0 0\n";
	print OUT "Pharmacophore feature\n";
	print OUT "AMBER ff99bsc0\n\n\n";
	print OUT "@<TRIPOS>ATOM\n";
	for(my$i =0; $i < scalar@{$X}; $i++)
	{
		my$ee = $i+1;
		if($$E[$i] <$S[1])			               #######Change the Cutoff accordingly
		{
			print OUT "      $ee ATM      $$X[$i]   $$Y[$i]   $$Z[$i] O.2       1 UNK    $$E[$i]\n";
		}
	}
	print OUT "@<TRIPOS>BOND\n";
	print OUT "@<TRIPOS>SUBSTRUCTURE\n";
    print OUT "1 UNK     1 RESIDUE           4 A     UNK     0 ROOT\n";
	
	close(OUT);
}


sub CliquePharm_Input
{
	my($X,$Y,$Z,$E,$Probe_Type,$name,$file) = @_;
	my@S = split(":",$Probe_Type);										## Split probe with MIF cutoff value 
	$name = lc($name);													## protein target name
	$Probe_Type = lc($S[0]);
	my@tmp1 = split("-",$file);											### Split input file name 
	my@tmp2 = split("_",$tmp1[2]);
	
	my$target_file = "energy_$name$tmp1[1]$tmp2[0]\_$S[0].pl";
	open(OUT,">../CliqueP_CathD/$target_file");                              #### Change Output Directory Name.
	
	for(my$i =0; $i < scalar@{$X}; $i++)
	{
		if($$E[$i] < $S[1])				        						#######Change the Cutoff accordingly
		{
			print OUT ("has_energy($name,$Probe_Type,$$X[$i],$$Y[$i],$$Z[$i],$$E[$i]).\n");
		}
	}
}
