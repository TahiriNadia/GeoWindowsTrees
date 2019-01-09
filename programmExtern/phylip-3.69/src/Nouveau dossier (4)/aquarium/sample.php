<?php
session_start();
include_once "conexion.php";

if (isset($_POST["user"])) {
    echo $_POST["user"];
}

// echo 'Current script owner: ' . get_current_user();

$date = date("d-m-Y");
$heure = date("H:i");
// exec("/home/tahiri_n/public_html/aquarium/file.pl & ") 

?>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<HTML>

<html xmlns="http://www.w3.org/1999/xhtml">
<!-- DW6 -->
<head>
<!-- Copyright 2005 Macromedia, Inc. All rights reserved. -->
<title>Aquarium - Bilan</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<link rel="stylesheet" href="mm_health_nutr.css" type="text/css" />
<script language="JavaScript" type="text/javascript">
//--------------- LOCALIZEABLE GLOBALS ---------------
var d=new Date();
var monthname=new Array("Janvier","Février","Mars","Avril","Mai","Juin","Juillet","Août","Septembre","Octobre","Novembre","Decembre");
//Ensure correct for language. English is "January 1, 2004"
var TODAY = d.getDate() + ", " + monthname[d.getMonth()] + " " + d.getFullYear();
// <input type="hidden" name="user2" value="<? echo "$user" ?>"> 
//---------------   END LOCALIZEABLE   ---------------


</script>

<?php
function fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe) {

	// 1 : on ouvre le fichier
	$monfichier = fopen('./'.$nomFile, 'w');

	//= affichage du contenu de la BD
	$sql = "SELECT ".$RegroupementPar.", COUNT(*) AS Frequentation FROM aquarium GROUP BY ".$RegroupementPar;
	$res = mysql_query($sql);

	$nb = mysql_numrows($res);  // on recupere le nombre d'enregistrements


	$sql2 = "SELECT ".$RegroupementPar.", COUNT(*) AS Frequentation FROM aquarium GROUP BY ".$RegroupementPar;
	$res2 = mysql_query($sql2);

	$nb2 = mysql_numrows($res2);  // on recupere le nombre d'enregistrements

	// 2 : on écrit la première ligne du fichier
	fputs($monfichier, "<chart>\n");
	fputs($monfichier, "\t<chart_data>\n");
	fputs($monfichier, "\t\t<row>\n");


	fputs($monfichier, "\t\t\t<null/>\n");

	$j = 0;
	while ($j < $nb2){ // parcours des resultats de la requete

		$ListeJourSemaine = mysql_result($res2, $j, $RegroupementPar);
		$Frequentation2 = mysql_result($res2, $j, "Frequentation");
		
		$contenu = $ListeJourSemaine."\n(".$Frequentation2.")";
		fputs($monfichier, "\t\t\t<string>$contenu</string>\n");
		$j++; 

	}

	fputs($monfichier, "\t\t</row>\n");

	fputs($monfichier, "\t\t<row>\n");
	fputs($monfichier, "\t\t\t<string>Fréquentation en fonction :".$TitreGraphe."</string>\n");
	$i = 0;
	while ($i < $nb){ // parcours des resultats de la requete

		$JourSemaine = mysql_result($res, $i, $RegroupementPar);
		$Frequentation = mysql_result($res, $i, "Frequentation");

		fputs($monfichier, "\t\t\t<number>$Frequentation</number>\n");
		
		$i++; 

	}

	fputs($monfichier, "\t\t</row>\n");
	fputs($monfichier, "\t</chart_data>\n");
	fputs($monfichier, "</chart>\n");
	  
	// 3 : quand on a fini de l'utiliser, on ferme le fichier
	fclose($monfichier);
}


function fichierDatePourGraphique($nomFile, $RegroupementPar, $TitreGraphe) {

	// 1 : on ouvre le fichier
	$monfichier = fopen('./'.$nomFile, 'w');

	//= affichage du contenu de la BD
	$sql = "SELECT ".$RegroupementPar.", COUNT(*) AS Frequentation FROM aquarium GROUP BY ".$RegroupementPar;
	$res = mysql_query($sql);

	$nb = mysql_numrows($res);  // on recupere le nombre d'enregistrements


	$sql2 = "SELECT ".$RegroupementPar.", COUNT(*) AS Frequentation FROM aquarium GROUP BY ".$RegroupementPar;
	$res2 = mysql_query($sql2);

	$nb2 = mysql_numrows($res2);  // on recupere le nombre d'enregistrements

	// 2 : on écrit la première ligne du fichier
	fputs($monfichier, "<chart>\n");
	fputs($monfichier, "\t<chart_data>\n");
	fputs($monfichier, "\t\t<row>\n");


	fputs($monfichier, "\t\t\t<null/>\n");
	
	// Création d'un tableau associatif des jours de la semaine.
	$TabJourSemaine = array('Lundi' => 0, 'Mardi' => 0, 'Mercredi' => 0, 'Jeudi' => 0, 'Vendredi' => 0, 'Samedi' => 0, 'Dimanche' => 0);
	
	$j = 0;
	while ($j < $nb2){ // parcours des resultats de la requete

		$ListeJourSemaine = mysql_result($res2, $j, $RegroupementPar);
		$Frequentation2 = mysql_result($res2, $j, "Frequentation");
		if($ListeJourSemaine=="Mon"){
			//$ListeJourSemaine="Lundi";
			$TabJourSemaine['Lundi'] = $Frequentation2;
		}else if($ListeJourSemaine=="Tue"){
			//$ListeJourSemaine="Mardi";
			$TabJourSemaine['Mardi'] = $Frequentation2;
		}else if($ListeJourSemaine=="Wed"){
			//$ListeJourSemaine="Mercredi";
			$TabJourSemaine['Mercredi'] = $Frequentation2;
		}else if($ListeJourSemaine=="Thu"){
			//$ListeJourSemaine="Jeudi";
			$TabJourSemaine['Jeudi'] = $Frequentation2;
		}else if($ListeJourSemaine=="Fri"){
			//$ListeJourSemaine="Vendredi";
			$TabJourSemaine['Vendredi'] = $Frequentation2;
		}else if($ListeJourSemaine=="Sat"){
			//$ListeJourSemaine="Samedi";
			$TabJourSemaine['Samedi'] = $Frequentation2;
		}else{
			//$ListeJourSemaine="Dimanche";
			$TabJourSemaine['Dimanche'] = $Frequentation2;
		}
		//$contenu = $ListeJourSemaine."\n(".$Frequentation2.")";
		//fputs($monfichier, "\t\t\t<string>$contenu</string>\n");
		$j++; 

	}
	
	foreach ($TabJourSemaine as $key => $value) {
		// Affichage frequentation par jour de la semaine
		if($value!=0){
			$contenu = $key."\n(".$value.")";
			fputs($monfichier, "\t\t\t<string>$contenu</string>\n");
		}
	}

	fputs($monfichier, "\t\t</row>\n");

	fputs($monfichier, "\t\t<row>\n");
	fputs($monfichier, "\t\t\t<string>Fréquentation en fonction :".$TitreGraphe."</string>\n");
	
	
	foreach ($TabJourSemaine as $key => $value) {
		// Affichage frequentation par jour de la semaine
		if($value!=0){
			fputs($monfichier, "\t\t\t<number>$value</number>\n");
		}
	}
	
	/*$i = 0;
	while ($i < $nb){ // parcours des resultats de la requete

		$JourSemaine = mysql_result($res, $i, $RegroupementPar);
		$Frequentation = mysql_result($res, $i, "Frequentation");

		fputs($monfichier, "\t\t\t<number>$Frequentation</number>\n");
		
		$i++; 

	}*/

	fputs($monfichier, "\t\t</row>\n");
	fputs($monfichier, "\t</chart_data>\n");
	fputs($monfichier, "</chart>\n");
	  
	// 3 : quand on a fini de l'utiliser, on ferme le fichier
	fclose($monfichier);
}
/*

$nomFile = "sample.xml";
$RegroupementPar = "JourSemaine";
fichierPourGraphique($nomFile, $RegroupementPar);

*/

?>


<script language="javascript">AC_FL_RunContent = 0;</script>
<script language="javascript"> DetectFlashVer = 0; </script>
<script src="AC_RunActiveContent.js" language="javascript"></script>
<script language="JavaScript" type="text/javascript">
<!--
var requiredMajorVersion = 10;
var requiredMinorVersion = 0;
var requiredRevision = 45;
-->
</script>
<BODY bgcolor="#FFFFFF">


</head>

<body bgcolor="#F4FFE4">
<table width="100%" border="0" cellspacing="0" cellpadding="0">
  <tr bgcolor="#D5EDB3">
    <td colspan="3" rowspan="2"><img src="header_bkgr_sos_template_teacher.jpg" width="751" height="139" /></td>
    <td height="50" colspan="3" id="logo" valign="bottom" align="center" nowrap="nowrap">Bienvenue &agrave; l'Aquarium </td>
    <td width="5">&nbsp;</td>
  </tr>

  <tr bgcolor="#D5EDB3">
    <td height="51" colspan="3" id="tagline" valign="top" align="center">D&eacute;partement d'Informatique UQAM </td>
	<td width="5">&nbsp;</td>
  </tr>

  <tr>
    <td colspan="7" bgcolor="#5C743D"><img src="mm_spacer.gif" alt="" width="1" height="2" border="0" /></td>
  </tr>

  <tr>
    <td colspan="7" bgcolor="#99CC66" background="mm_dashed_line.gif"><img src="mm_dashed_line.gif" alt="line decor" width="4" height="3" border="0" /></td>
  </tr>

  <tr bgcolor="#99CC66">
  	<td colspan="7" id="dateformat" height="20">&nbsp;&nbsp;<script language="JavaScript" type="text/javascript">
      document.write(TODAY);	</script>	</td>
  </tr>
  <tr>
    <td colspan="7" bgcolor="#99CC66" background="mm_dashed_line.gif"><img src="mm_dashed_line.gif" alt="line decor" width="4" height="3" border="0" /></td>
  </tr>

  <tr>
    <td colspan="7" bgcolor="#5C743D"><img src="mm_spacer.gif" alt="" width="1" height="2" border="0" /></td>
  </tr>

 <tr>
    <td width="169" valign="top" bgcolor="#5C743D">
	<table border="0" cellspacing="0" cellpadding="0" width="165" id="navigation">
        <tr>
          <td width="165">&nbsp;<br />
		 &nbsp;<br /></td>
        </tr>
        <tr>
          <td width="165"><a href="index.html">Page d'accueil </a><a href="javascript:;" class="navText"></a></td>
        </tr>
        <tr>
          <td width="165"><a href="employes.php">Employ&eacute;s </a><a href="javascript:;" class="navText"></a></td>
        </tr>
		<tr>
          <td width="165"><a href="logout.php">Deconnection </a><a href="javascript:;" class="navText"></a></td>
        </tr>
        <tr>
          <td width="165"><a href="calendar.html">Calendrier</a><a href="javascript:;" class="navText"></a></td>
        </tr>
        <tr>
          <td width="165"><a href="contact.html">Contact</a><a href="javascript:;" class="navText"></a></td>
        </tr>
      </table>
 	 <br />
  	&nbsp;<br />
  	&nbsp;<br />
  	&nbsp;<br /> 	</td>
    <td width="50"><img src="mm_spacer.gif" alt="" width="50" height="1" border="0" /></td>
    <td colspan="2" valign="top"><img src="mm_spacer.gif" alt="" width="305" height="1" border="0" /><br />
	&nbsp;<br />
	
<?php
$nomFile = "JourSemaine.xml";
$TitreGraphe = " du jour de la semaine.";
$RegroupementPar = "JourSemaine";
fichierDatePourGraphique($nomFile, $RegroupementPar, $TitreGraphe);


$nomFile = "PeriodeJournee.xml";
$TitreGraphe = " de la période de la journée.";
$RegroupementPar = "PeriodeJournee";
fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe);


$nomFile = "CodeProgramme.xml";
$TitreGraphe = " du code de programme.";
$RegroupementPar = "CodeProgramme";
fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe);

$nomFile = "CoursSuivi.xml";
$TitreGraphe = " du cours suivi.";
$RegroupementPar = "CoursSuivi";
fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe);

$nomFile = "NomEmploye.xml";
$TitreGraphe = " du nom de l'employé(e).";
$RegroupementPar = "NomEmploye";
fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe);

$nomFile = "NumeroSemaine.xml";
$TitreGraphe = " du numéro de la semaine.";
$RegroupementPar = "DATE_FORMAT(Date,'%u')";
fichierPourGraphique($nomFile, $RegroupementPar, $TitreGraphe);




//setlocale (LC_TIME, 'fr_FR.utf8','fra'); 
//echo (strftime("%A %d %B")); 

//echo "Session Automne 2013: La fréquentation totale de l'aquarium datant du ."(strftime("%A %d %B"))." est de : ".$res;


?>

<script language="JavaScript" type="text/javascript">
	VoirGraphe('JourSemaine.xml');
	VoirGraphe('PeriodeJournee.xml');
	VoirGraphe('CodeProgramme.xml');
	VoirGraphe('CoursSuivi.xml');
	VoirGraphe('NomEmploye.xml');
	VoirGraphe('NumeroSemaine.xml');
	
<!--
function VoirGraphe(nomFile) {
	this.nomFile = nomFile;

if (AC_FL_RunContent == 0 || DetectFlashVer == 0) {
	alert("This page requires AC_RunActiveContent.js.");
} else {
	var hasRightVersion = DetectFlashVer(requiredMajorVersion, requiredMinorVersion, requiredRevision);
	if(hasRightVersion) { 
		AC_FL_RunContent(
			'codebase', ' ',
			'width', '400',
			'height', '250',
			'scale', 'noscale',
			'salign', 'TL',
			'bgcolor', '#777788',
			'wmode', 'opaque',
			'movie', 'charts',
			'src', 'charts',
			'FlashVars', 'library_path=charts_library&xml_source='+this.nomFile, 
			'id', 'my_chart',
			'name', 'my_chart',
			'menu', 'true',
			'allowFullScreen', 'true',
			'allowScriptAccess','sameDomain',
			'quality', 'high',
			'align', 'middle',
			' ', ' ',
			'play', 'true',
			'devicefont', 'false'
			); 
	} else { 
		var alternateContent = 'This content requires the Adobe Flash Player. '
		+ '<u><a href=http://www.macromedia.com/go/getflash/>Get Flash</a></u>.';
		document.write(alternateContent); 
	}
}
}
// -->
</script>
<noscript>
	<P>This content requires JavaScript.</P>
</noscript>




&nbsp;<br />
	<table border="0" cellspacing="0" cellpadding="0" width="502">
	<tr>
		<td>
			<?php 

				//= affichage du contenu de la BD
				$sql = "SELECT COUNT(*) AS Frequentation_Totale FROM aquarium WHERE Session LIKE '%2014%'";
				$res ="";
				$res = mysql_query($sql);
				$nb=0;
				$nb = mysql_numrows($res);  // on recupere le nombre d'enregistrements
				$Freq=0;
				
				$j = 0;
				while ($j < $nb){ // parcours des resultats de la requete

					$Freq = mysql_result($res, $j, "Frequentation_Totale");
					/*var_dump($Freq);
					$Freq = trim($Freq);
					$Freq = rtrim($Freq); 
					var_dump($Freq);*/
					$j++; 

				}
				
				
				$jour = array("Dimanche","Lundi","Mardi","Mercredi","Jeudi","Vendredi","Samedi"); 
				$mois = array("","Janvier","Février","Mars","Avril","Mai","Juin","Juillet","Août","Septembre","Octobre","Novembre","Décembre"); 
				$datefr = $jour[date("w")]." ".date("d")." ".$mois[date("n")]." ".date("Y"); 

				echo "<BR><u>Session Hiver 2014</u> : <BR>La fréquentation totale de l'aquarium datant du ".$datefr." est de : ".$Freq;
				//echo "Nous sommes le ". $datefr; 
			?> 
		</td>
		
	</tr>
	
	
		
      </table>
	 <br />
	&nbsp;<br />	</td>
    <td width="6">&nbsp;</td>
        <td width="4" valign="top">&nbsp;</td>
	<td width="5">&nbsp;</td>
  </tr>
  
  <tr>
    <td width="169">&nbsp;</td>
    <td width="50">&nbsp;</td>
    <td width="534">&nbsp;</td>
    <td width="256">&nbsp;</td>
    <td width="6">&nbsp;</td>
    <td width="4">&nbsp;</td>
	<td width="5">&nbsp;</td>
  </tr>
</table>

</BODY>
</HTML>
