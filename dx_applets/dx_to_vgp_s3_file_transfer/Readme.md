# DNAnexus to VGP S3 Exporter

## What Does this App Do?

This app will upload file(s) to the official VGP S3 bucket.

## How Does this App Work?

This app uses AWS CLI (Amazon Web Services Command Line Interface) to upload files to the VGP's S3 bucket. 
Once the files are transferred to the S3 bucket, the original file is removed from the DNAnexus platform and
is replaced with a link to the new file in the S3 bucket.

The files will be organized in the S3 bucket by species name.  By default, this app will look for a property
named `species_name` on the project the job is running in.  If the app finds this property, then it will 
use this as the species name and transfer the given files to this species' folder in the S3 bucket.  If no
`species_name` property is found on this project, then it will look at the value provided as an input parameter
when the app was run.  The app will fail if it can not find a species name in one of these locations.
