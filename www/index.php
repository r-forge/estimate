<!-- $Id$ -->
<!-- This is the project specific website -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

?>

<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title><?php echo $group_name; ?></title>
  <link rel="stylesheet"
        type="text/css"
        href="http://<?php echo $themeroot; ?>styles/estilo1.css" />
</head>
<body>
  <header>
    <!-- R-Forge Logo -->
    <a href="http://r-forge.r-project.org/">
      <img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png"
           border="0"
           alt="R-Forge Logo" />
    </a>
  </header>

  <div class="content">
    <p>
      <abbr title="Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data">ESTIMATE</abbr> is a tool for predicting tumor
      purity, and the presence of infiltrating stromal/immune cells in tumor
      tissues using gene expression data. Its algorithm is based on single
      sample Gene Set Enrichment Analysis, and generates three scores:
      <ol>
        <li>stromal score (that captures the presence of stroma in tumor tissue)</li>
        <li>immune score (that represents the infiltration of immune cells in tumor tissue)</li>
        <li>estimate score (that infers tumor purity)</li>
      </ol>
    </p>

    <p>
      Click
      <a href="http://bioinformatics.mdanderson.org/main/ESTIMATE:Overview">
      here</a> for the software's official
      <strong><abbr title="MD Anderson Cancer Center">MDACC</abbr> jump
      page</strong>, which includes a link to the ESTIMATE scores for all
      <abbr title="The Cancer Genome Atlas">TCGA</abbr> tumor types.
      <br />

      Click
      <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a> for the R-Forge <strong>project summary page</strong>.
    </p>
  </div>
</body>
</html>

