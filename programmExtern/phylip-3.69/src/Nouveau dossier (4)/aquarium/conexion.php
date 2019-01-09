<?php
// datos para la coneccion a mysql
define('DB_SERVER','wozniak.uqam.ca');
define('DB_NAME','aquarium');
define('DB_USER','tahiri_n');
define('DB_PASS','123456');

$con = mysql_connect(DB_SERVER,DB_USER,DB_PASS);
mysql_select_db(DB_NAME,$con);
?>

