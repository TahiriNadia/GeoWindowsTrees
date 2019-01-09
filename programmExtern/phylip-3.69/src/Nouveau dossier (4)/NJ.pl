#!/bin/perl
#=========================================================================================
#= Apllel du programme Neigbor joining
#= Auteur : Tahiri Nadia
#= Date	  : Automne 2014
#=========================================================================================

#== Les inclusions
use strict;
use warnings;

#== Appel de NJ
my $cmd = "./neighbor <parametres.txt";
system($cmd);


my $nom_fichier = "outtree";

#== ouverture du fichier
open (IN , "$nom_fichier") || die ("$!");
open (OUT , ">>out_trees") || die ("$!");


#== lecture du contenu du fichier
while (my $ligne = <IN>){
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


#== Message de fin de script
print STDOUT "\n\nFin normale du script $0\n";