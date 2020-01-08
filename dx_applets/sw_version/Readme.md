# Cloud Workstation App

## What does this app do?
This app sets up a cloud workstation which you can access by running the app with the **--ssh** or **--allow-ssh** flags.

## What are typical use cases for this app?
This app can be used as a workstation inside of the DNAnexus cloud platform.  By running the app with **--ssh** or
**--allow-ssh**, users can login to a machine inside of the DNAnexus cloud platform.  From there, users can upload/download
data to/from the project in which the app is run, perform data analysis, and install additional packages from sources such as 
apt, cran, pip, github, etc.

One note: in order to access files stored in the project in which this app is being run, users must provide the project-identifier.
This can easily be done with a special environment variable called $DX\_PROJECT\_CONTEXT\_ID.  For instance, to download a file called
myreads.fastq.gz from the parent project, users would simply run 
```dx download $DX_PROJECT_CONTEXT_ID:/myreads.fastq.gz```

## What are the inputs?
The user can provide a maximum session length value.  After this amount of time has passed, the workstation will automatically 
shut-down.  Timeout is provided using suffixes s, m, h, d, w, M, y.

During a session, users can check how much time remains until the session times out by running
```dx-get-timeout```
Users can reset the timeout by running
```dx-set-timeout <timeout>```
using the suffixes s, m, h, d, w, M, y.

Additionally, users can provide a list of files to download at the app startup.  These files will be downloaded to the home
directory automatically, allowing easy access for data analysis.

Finally, users can create snapshots of their cloud environment by running
```dx-snapshot <snapshot name (optional)>```
or
```dx-close-and-save-notebook <snapshot name (optional)>```
These will generate .snapshot files which can be provided as input to future cloud-workstation runs
to pick-up where you left off.

