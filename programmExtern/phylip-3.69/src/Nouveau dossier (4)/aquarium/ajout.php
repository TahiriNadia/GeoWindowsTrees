<?php
session_start();
include_once "conexion.php";

if (isset($_POST["user"])) {
    echo $_POST["user"];
}

// echo 'Current script owner: ' . get_current_user();

$date = date("d-m-Y");
$heure = date("H:i");
?>




<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<!-- DW6 -->
<head>
<!-- Copyright 2005 Macromedia, Inc. All rights reserved. -->
<title>Aquarium - Ajout</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<link rel="stylesheet" href="mm_health_nutr.css" type="text/css" />
<script language="JavaScript" type="text/javascript">
//--------------- LOCALIZEABLE GLOBALS ---------------
var d=new Date();
var monthname=new Array("January","February","March","April","May","June","July","August","September","October","November","December");
//Ensure correct for language. English is "January 1, 2004"
var TODAY = monthname[d.getMonth()] + " " + d.getDate() + ", " + d.getFullYear();
// <input type="hidden" name="user2" value="<? echo "$user" ?>"> 
//---------------   END LOCALIZEABLE   ---------------
</script>
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
  	<td colspan="7" id="dateformat" height="20">&nbsp;&nbsp;<?php echo date("F d, Y"); ?>	</td>
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
    <td colspan="2" valign="top"><img src="mm_spacer.gif" alt="" width="305" height="1" border="0" /><br /><br />
	
	
	<form action="" method="post" class="login">
		<!--<div><p id="periode"><label>P&eacute;riode de la journ&eacute;e : </label><p align="right">
			<input name="periode" type="radio" value="AM" checked="checked" / >
			AM
			<input type="radio" name="periode" value="PM" />PM
			<input type="radio" name="periode" value="SOIR" />SOIR
		</p>
		</div>-->
		<div><label>Code programme : </label><input name="programme" type="text" ></div>
		<div><label>Cours suivi : </label><input name="cours" type="text" ></div>
		<div><label>Dur&eacute;e de l'aide (en minutes) : </label><input name="duree" type="text" ></div>		
		<div><label>Type probl&egrave;me : </label><br><br><br><textarea id="pb" name="probleme" row="200" cols="55"></textarea></div>
		<div><label>Commentaires : </label><br><br><br><textarea id="comm" name="commentaires" row="200" cols="55" ></textarea></div>
		<div><input name="ajouter" type="submit" value="ajouter"><?php $valide=1; ?></div>
	</form>
	
	<?php
			if(isset($_POST['duree']) && isset($_POST['programme']) && isset($_POST['cours']) && isset($_POST['probleme']) ){
				//= insertion des donnees
				// $nom = $_POST['nom'];
				$nom = $_SESSION['nom'];
				//$periode = $_POST['periode'];
				$duree = $_POST['duree'];
				$programme = $_POST['programme'];
				$cours = $_POST['cours'];
				$probleme = $_POST['probleme'];
				$commentaires = $_POST['commentaires'];
				
				$periode = "";
			}
			
			// if(isset($_POST['nom']) && isset($_POST['periode']) && isset($_POST['duree']) && isset($_POST['programme']) && isset($_POST['cours']) && isset($_POST['probleme']) && isset($_POST['commentaires']) ){
			if(!empty ($duree) && !empty ($programme) && !empty ($cours) && !empty ($probleme) ){
				// echo $_POST['nom'];	
				// echo $_POST['periode'];
				// echo $_POST['duree'];
				// echo $_POST['programme'];
				// echo $_POST['cours'];
				// echo $_POST['probleme'];
				// echo $_POST['commentaires'];  
				$sql = "INSERT INTO  `aquarium` (`NomEmploye`, `Date`, `JourSemaine`, `PeriodeJournee`, `DureeAide`, `CodeProgramme`, `CoursSuivi`, `TypePb`, `Commentaire`) VALUES('" . $nom . "', CURDATE(), DATE_FORMAT(CURDATE(),'%a'), '" . $periode . "', " . $duree . ", '".$programme."','" . $cours . "','" . $probleme ."','" . $commentaires ."');";
				
				// echo $sql;
				$res = mysql_query($sql);
				
				$nom = "";
				$periode = "";
				$duree = "";
				$programme = "";
				$cours = "";
				$probleme = "";
				$commentaires = ""; 
				
				// unset($_POST['nom']);
				// unset($_POST['periode']);
				// unset($_POST['duree']);
				// unset($_POST['programme']);
				// unset($_POST['cours']);
				// unset($_POST['probleme']);
				// unset($_POST['commentaires']); 
				
				echo "<font align=right color=red>Enregistrement effectu&eacute;!</font>";
			}else{
				if(isset($_POST['ajouter'])){
					echo "<font align=right color=red>Veuillez compl&eacute;ter les champs suivants : dur&eacute;e, programme, cours et probl&egrave;me!</font>";
				}else{
				}
				
			}
			$valide=0;
	?>
	<style type="text/css">
*{
    font-size: 14px;
}
body{
background:#aaa;
}
form.login {
    background: none repeat scroll 0 0 #F1F1F1;
    border: 1px solid #DDDDDD;
    font-family: sans-serif;
    margin: 0 auto;
    padding: 20px;
    width: 500px;
    box-shadow:0px 0px 20px black;
    border-radius:10px;
}
form.login div {
    margin-bottom: 15px;
    overflow: hidden;
}
form.login div label {
    display: block;
    float: left;
    line-height: 25px;
}
form.login div input[type="text"], form.login div input[type="password"] {
    border: 1px solid #DCDCDC;
    float: right;
    padding: 4px;
}
form.login div input[type="submit"] {
    background: none repeat scroll 0 0 #DEDEDE;
    border: 1px solid #C6C6C6;
    float: right;
    font-weight: bold;
    padding: 4px 20px;
}
.error{
    color: red;
    font-weight: bold;
    margin: 10px;
    text-align: center;
}
</style>
	&nbsp;<br />
	<table border="0" cellspacing="0" cellpadding="0" width="502">
        
		

		
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
</body>
</html>
