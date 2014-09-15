##### Data Help #####
Data="<p>The Data Import module assists in bringing data into the application for downstream analysis.  Two data files are required for analysis: an InCell derived file, and a REMP Platemap file.  For convenience, these files can be uploaded directly to the app using the file upload functions; however, there is a maximum upload size of 45 Mb.  To work around that limitation, files can be placed in the Shiny Dropbox folder on the storage server and accessed from within the app using the dropdowns.</P>
<hr>
<h3>Inputs</h3>
<ul>
  <li><p><strong>Import InCell File</strong> Choose an InCell file to analyze using either the file upload tool or, if the file is already in the Shiny Dropbox, choose the file from the dropdown.  There is a 45 Mb maximum file size for uploading files.  InCell files must be saved in CSV format.  The file parser will find the header rows for each data level and separate them into distinct tables.</p></li>
  <li><p><strong>Import REMP PlateLookup File</strong> To annotate the data with the compounds and concentrations, a REMP PlateLookup file is required.  This file should contain only data for the compound and control plates that were used to treat the assay plate in the InCell file.  If mutliple compound or control plates are detected in the file, a reasonable attempt will be made to include annotations from one compound plate and one control plate only.</p></li>
  <li><p><strong>Check for New Files</strong> If the file you are looking for is not listed in one of the dropdown boxes, and it was recently placed in the Shiny Dropbox, clicking this button will rescan the folder and add the file to the proper dropdown.</p></li>
  <br/>
  <li><p><strong>Control Columns</strong> The user must specify which columns received compounds from the control plate.  Control columns can be specified by listing each column sparated by a comma (ie 21,33,23,24) or listing the minimum and maximum columns separated by a hyphen (ie 21-24).  The default is to use columns 21-24 as control columns.</p></li>
  <li><p><strong>Set Positive Control</strong> From the list of compounds in the REMP PlateLookup file, the user must choose a compound to use as the positive control in later analysis.  In addition to choosing a positive control compound, the user must also choose the concentration of that compound to use as the positive control.  The default concentration is the maximal concentration for that compound.</p></li>
</ul>
<hr>
<h3>Outputs</h3>
<ul>
    <li><p><strong>Summary Statistics</strong> Low-level statistics about the plate derived from the data files loaded.</p>
      <ul>
          <li><em>Data Files:</em> The file names of the InCell and REMP PlateLookup files loaded</li>
          <li><em>Total Wells:</em> The number of wells of data loaded from the InCell file</li>
          <li><em>Total Cells Counted:</em> The total number of cell for which there is cell-level data</li>
          <li><em>Total Fields:</em> The total number of fields loaded</li>
          <li><em>Mean Fields Per Well:</em> The average number of fields per well.  This number is calculated from the total number of fields and the number of wells.</li>
        </ul>
      </li>
      <li><strong>Data File Download</strong> Each data table is available to download individually in CSV format.  Simply click the download button for the data table needed.  Data table headers are compressed to one row by adding a hyphen between lines.</li>
</ul>
<br/>"



##### QC Help #####
QC="<p>The Quality Control module is designed so that glaring problems in the dataset can be quickly and easily identified.  A large heatmap displays the values of any feature for quick visual inspection, while a small histogram shows the distribution of the feature across the plate.  Click on any well of the heatmap to generate a boxplot showing the distribution the selected feature within that well.</P>
<hr>
<h3>Inputs</h3>
<ul>
  <li><p><strong>Select a Feature</strong> This dropdown input is populated from the headings found in the InCell file that is loaded.  When 'Cell Number' is selected, the intrawell distribution boxplot will not be drawn.</p></li>
</ul>
<hr>
<h3>Outputs</h3>
<ul>
  <li><p><strong>Summary Statistics</strong> Descriptive statistics of the entire plate:</p>
      <ul>
        <li><em>Min:</em> The lowest value in the dataset</li>
        <li><em>1st Quartile:</em> The data value marking the boundary between the lowest 25% of values in the dataset and the rest of the data</li>
        <li><em>Median:</em> The data value that separates the highest half of the data from the lowest half</li>
        <li><em>Mean:</em> The central value of the dataset</li>
        <li><em>3rd Quartile:</em> The data value marking the boundary between the lowest 75% of values in the dataset and the rest of the data</li>
        <li><em>Max:</em> The highest value in the dataset</li>
      </ul>
    </li>
    <li><p><strong>Heatmap</strong> A color-coded representation of the values of each well of the dataset arrayed in the format of a 384-well plate.  The X and Y  axes represent physical coordinates of the corresponding well in the plate.  The color scale represents the range of values for the chosen feature. Hover the cursor over a well to view it's coordinates and the value of the chosen feature averaged across all cells in the well.  Clicking the well produces a new boxplot output.</p></li>
    <li><p><strong>Histogram</strong> The distribution of the chosen feature across the entire plate.  Useful for quickly determining whether data is skewed or otherwise not normally distributed</p></li>
    <li><p><strong>Boxplot</strong> A box-and-whiskers plot displaying the distribution of the selected feature within the selected well.  Click on a well in the heatmap to generate or update the boxplot.  Each imaging field is included as a seprate category in the plot; the DMSO treated control well distribution is appended to the fields of the selected well as an additional category.</p></li>
</ul>
<br/>"


##### Any|Any Help #####
AA="<p>The Any|Any module facilitates plotting any feature versus any other feature.  Feature values can be log transformed and/or thresholded.  The plot can be annotated with a single horizontal and a single vertical line.  Lines can be labeled and placed at any value on the axis.  Smaller plots, below the main feature plot, assist in setting thresholds by showing the distribution of positive and negative controls.</P>
<hr>
<h3>Inputs</h3>
<ul>
  <li><p><strong>Y/X Axis Feature</strong> The feature to plot on the specified axis.  The dropdown is populated from the column heading in the InCell file provided.</p></li>
  <li><p><strong>Y/X Axis Transform</strong> Which type of transformation to apply to the data.</p></li>
  <li><p><strong>Y/X Axis Threshold</strong> Whether to apply a threshold to the data.  If an option other than 'None' is selected the values plotted will be the percentage above or below the value set with the slider.</p></li>
  <li><p><strong>Y/X Axis Threshold Slider</strong> The threshold value selector.  Use the Plus and Minus buttons to adjust the slider value in increments of one; drag the slider for coarse control.</p></li>
  <li><p><strong>Add Horizontal/Vertical Line</strong> Whether to show a horizontal or vertical line on the plot.</p></li>
  <li><p><strong>Label</strong> The text to use as a label for the specified line.</p></li>
  <li><p><strong>Value</strong> The value at which to draw the line.</p></li>
</ul>
<hr>
<h3>Outputs</h3>
<ul>
    <li><p><strong>Any|Any Plot</strong> A scatter plot showing the value of each well of the plate in the coordiantes of the features chosen.  Each compound is plotted as a series whith the symbol size corresponding to concentration - larger symbols means higher concentration.  Series can be hidden by clicking the name of the compound in the legend.  If the legend length exceeds the height of the plot the legend will be scrollable.  Hovering the cursor over a symbol on the plot will highlight that point and the rest of it's series.</p></li>
    <li><p><strong>Y/X Axis Thresholding Tool</strong> The positive and negative control values of the feature chosen for the specified axis are plotted to show the distribution of the controls and assist in choosing a threshold value that can differentiate them.  When an option other than 'None' is chosen for the Axis Threshold this plot will show either the percent above or the percent below the Threshold Slider value.  An optimal threshold value will result in all negative controls at 0% and all positive controls at 100%.</p></li>
</ul>
<br/>"


##### Any|Conc Help #####
AC="<p>The Any|Conc module facilitates plotting any feature versus the compound concentration.  Feature values can be log transformed and/or thresholded.  A single, smaller plot, below the main feature plot, assist in setting thresholds by showing the distribution of positive and negative controls.</P>
<hr>
<h3>Inputs</h3>
<ul>
  <li><p><strong>Y Axis Feature</strong> The feature to plot on the specified axis.  The dropdown is populated from the column heading in the InCell file provided.</p></li>
  <li><p><strong>Y Axis Transform</strong> Which type of transformation to apply to the data.</p></li>
  <li><p><strong>Y Axis Threshold</strong> Whether to apply a threshold to the data.  If an option other than 'None' is selected the values plotted will be the percentage above or below the value set with the slider.</p></li>
  <li><p><strong>Y Axis Threshold Slider</strong> The threshold value selector.  Use the Plus and Minus buttons to adjust the slider value in increments of one; drag the slider for coarse control.</p></li>
</ul>
<hr>
<h3>Outputs</h3>
<ul>
    <li><p><strong>Any|Conc Plot</strong> A scatter plot showing the value of each well of the plate in the coordiantes of the feature chosen over the concentration of compound.  Each compound is plotted as a distinct series.  Series can be hidden by clicking the name of the compound in the legend.  If the legend length exceeds the height of the plot the legend will be scrollable.  Hovering the cursor over a symbol on the plot will highlight that point and the rest of it's series.</p></li>
    <li><p><strong>Y Axis Thresholding Tool</strong> The positive and negative control values of the feature chosen are plotted to show the distribution of the controls and assist in choosing a threshold value that can differentiate them.  When an option other than 'None' is chosen for the Axis Threshold this plot will show either the percent above or the percent below the Threshold Slider value.  An optimal threshold value will result in all negative controls at 0% and all positive controls at 100%.</p></li>
</ul>
<br/>"


##### Feature Help #####
Feat="<p>The Feature module is designed to give an overview of the activity of a compound on a feature of interest.  By setting a feature and scrolling through the compounds, the user can quickly identify the compounds that are affecting the chosen feature.  A large, main plot shows the distribution of feature values for each concentration of the compound tested.  Shaded series show how the current compound compares to positive and negative controls.  For each concentration, statistical tests are performed to quantify any difference from the negative control.  These values are z-score transformed and plotted in a smaller plot below the main plot.  From this analysis, the area under each statistic's curve is calculated to give an overall difference from the negative control.</P>
<hr>
<h3>Inputs</h3>
<ul>
  <li><p><strong>Select a Feature</strong> The feature for which to plot distributions. This dropdown input is populated from the headings found in the InCell file that is loaded.  When 'Cell Number' is selected, no plots will be drawn.</p></li>
  <li><p><strong>Select a Compound</strong> The compound for which to plot distributions.  The values in the dropdown are populated from the REMP PlateLookup file that was selected in the 'Data' tab.  Use the previous and next buttons below the dropdown to advance or rollback the compound shown o nthe plot.</p></li>
  <li><p><strong>Log Transform Feature Values</strong> Which type of transformation to apply to the data.  Some features are better viewed on a logarithmic axis.</p></li>
</ul>
<hr>
<h3>Outputs</h3>
<ul>
    <li><p><strong>Distribution Plot</strong> Line plot showing the kernel density estimations for each concentration of the chosen compound.  The shaded series show kernel density estimations for negative and positive controls.  By default, every other series is hidden to make the plot less complicated.  To show the series, simply click the series name in the legend.  To show a zoomed view of the plot, click and drag to select an area of the plot to zoom in on.  In zoom view, a rest button will appeart to reset the plot.</p></li>
    <li><p><strong>Statistic Plot</strong> Comparisons between each concentration and the negative control are tested for statistical significance.  The test statistics generated are z-score normalized using the median and median absolute deviation of the negative controls.  The normalized test statistics are plotted over concentration to produce this plot.  Included as a series on this plot is the z-score normalized cell number for each well.  This enabled changes in the number of cells to be quickly identified as a confounder in the analysis.</p></li>
    <li><p><strong>Area Under the Cure Plot</strong> The shaded area under the statistic curves on the Statistic Plot are calculated and represented by bars in this plot.  Hover over the bar to view it's exact value.</p></li>
</ul>
<br/>
"


##### Action Help #####
Act=""
