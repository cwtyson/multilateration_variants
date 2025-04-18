# Tracking small animals in complex landscapes: a comparison of localisation workflows for automated radio telemetry systems

[https://doi.org/10.5061/dryad.3ffbg79t3](https://doi.org/10.5061/dryad.3ffbg79t3)

## Description of the data and file structure

This dataset includes all data and necessary code in R scripts to produce the results of the associated paper.

### Files and variables

#### File: tags\_id.csv

**Description:** Table with information of tags deployment on hummingbirds.

##### Variables

* tag_id: Unique identifier for tags
* tag_type: Type of tag, LifeTag or PowerTag
* tag_attachment_method: Tag attachment method, harness or glue-on
* attached: Whether tag has been deployed on an animal or not
* date_m_d_y: Date of tag deployment, in month, day and year
* species: Latin name of species tag was deployed on
* id: Unique identifier for individual in color code
* hour: Time of tag deployment

#### File: Nodes.csv

**Description:** Table with information on node location.

##### Variables

* Own_NodeId: Own unique node identifier
* NodeId: Manufacturer's unique node identifier
* Latitude: Latitude in decimal degrees
* Longitude: Longitude in decimal degrees
* NodeUTMx: UTM longitudinal location, projected EPSG
* NodeUTMy: UTM latitudinal location, projected EPSG

#### File: time\_diff\_matched\_df.csv

**Description:** Table with matched signal strength values for tags used in calibration trials.

##### Variables

* TagId: Unique identifier for tags
* Time.local: Time in GMT
* NodeId: Manufacturer's unique node identifier
* TagRSSI: Signal strength measured in decibels
* time_format: Formatted time to switch to local time
* node_local_time: Local time
* unique_trial_id: Unique identifier for calibration trial
* lon: Longitude in decimal degrees
* lat: Latitude in decimal degrees
* elevation: Altitude in metres above sea level
* test_height: Height of calibration trial, either ground (low), mid or drone (high)
* test_local_time: Local time of start of calibration trial
* tag_type: Type of tag, LifeTag or PowerTag
* start_time: Start time of calibration trial
* end_time: End time of calibration trial
* reloc_id: Unique identifier to location estimate
* time_difference: Difference in time between beep data from radio signal and measured time in calibration trial (measured in seconds)
* abs_time_difference: Correction of time difference for absolute values (measured in seconds)

#### File: beep\_data.zip

**Description:** Compressed folder with radio signal strength reads during the time that hummingbirds flew over the grid.

## Code/software

All analyses were carried out using R. The following files for R scripts are included:

1.decay_function.R: Calculate exponential decay function from calibration data

2.nls_multilateration.R: Estimate locations through different multilateration options

3.error_analyses.R: Calculate error in metres of estimates to known locations

4.hummer*_*rss_to_tracks.R: Convert signal reads to movement tracks for tracked hummingbirds

5.Second_smoothing_after_multilat.R: Smoothing options after multilateration

6.Smoothing_rss_before_multilat.R: Smoothing options before multilateration

## Access information

Other publicly accessible locations of the data:

* Data on ground elevation was obtained from the GLO-30 Copernicus digital elevation model (30 m resolution) of the European Space Agency (ESA, 2024).
