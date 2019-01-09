#!/bin/perl
#=========================================================================================
#= Apllel du programme Neigbor joining
#= Auteur : Tahiri Nadia
#= Date	  : Automne 2014
#=========================================================================================

#== Les inclusions
use strict;
use warnings;

#== découpement des matrices sur plusieurs fichiers
open (INPUT , "matrices.txt") || die ("$!");

my $ligne;
while($ligne = <INPUT>){
	chomp($ligne);
	open (OUTPUT , ">mat.txt") || die ("$!");
	
	if($ligne=~/\d/){
		my $nb = $ligne + 0;
		print STDOUT "$ligne\n";
		print OUTPUT "$ligne\n";
		for(my $i=1;$i<=$nb;$i++){
			$ligne = <INPUT>;
			chomp($ligne);
			print STDOUT "$ligne\n";
			print OUTPUT "$ligne\n";
		}
	}
	close(OUTPUT);
	#== Appel de NJ
	my $cmd = "./neighbor <parametres.txt";
	system($cmd);
	
	#== ouverture du fichier
	open (IN , "outtree") || die ("$!");
	open (OUT , ">>out_trees") || die ("$!");
	
	#== lecture du contenu du fichier
	while ($ligne = <IN>){
		chomp($ligne);
		print STDOUT "$ligne";
		print OUT "$ligne";
	}

	print STDOUT "\n";
	print OUT "\n";

	#== fermeture du fichier
	close IN;

	$cmd = "rm outfile";
	system($cmd);

	$cmd = "rm outtree";
	system($cmd);
}

#== Fermeture des fichiers
close(IN);
close(OUT);
close(INPUT);

#== Message de fin de script
print STDOUT "\n\nFin normale du script $0\n";