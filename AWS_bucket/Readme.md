# VGP Upload Script

## Dependencies

* boto3 (pip install boto3)
* dx-toolkit (https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK)
* aws cli (https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html)
* python2.7

Please authenticate your user using dx toolkit before running this script, by running:
`dx login` or `dx login --token [yourauthtoken]`

The script takes as input an AWS profile name. Please configure the profile using the AWS CLI instructions: https://docs.aws.amazon.com/cli/latest/userguide/cli-multiple-profiles.html

The script assumes that the DNAnexus drive and AWS profile will have the same name, so the use of named profiles is recommended.

The bucket name is hardcoded into the script under "VGP_BUCKET"

## Running

The script is run using python:
`python vgp_upload_to_cloud.py /path/to/file --species-name [species_name] --species-id [species_id] --datatype [datatype] --profile [profile]`

Allowed datatypes are:
* "10x" for 10X Genomics
* "pacbio" for Pacific Biosciences
* "arima", "phase" or "hic" for Hi-C (depending on data source)
* "bionano" for Bionano Genomics

The script will create the following directory structure based on provided inputs:
```
.
└── species
    └── <species_name>
       └── <species_id>
            ├── genomic_data
            │   ├── bionano
            │   │   ├── testfile1.txt
            │   │   ├── testfile2.txt
            │   │   ├── testfile3.txt
            │   │   └── testfile4.txt
            │   └── pacbio
            │       ├── testfile1.txt
            │       ├── testfile2.txt
            │       ├── testfile3.txt
            │       └── testfile4.txt
            └── transcriptomic_data
                └── <tissue_id>
                    ├── illumina
                    │   ├── testfile1.txt
                    │   ├── testfile2.txt
                    │   ├── testfile3.txt
                    │   └── testfile4.txt
                    └── pacbio
                        ├── testfile1.txt
                        ├── testfile2.txt
                        ├── testfile3.txt
                        └── testfile4.txt
```

In DNAnexus, the directory structure is replicated with 'species_id' serving as the Project name. The 'species_name' is not used but is added as a tag to all files.

## DX Drives and Symlinks

The script will automatically create a DNAnexus drive object with the specified AWS profile name and symlinks to associate file-ids to AWS S3 files. A drive is required to create symlinks. Below are API calls to do this manually.

Adding an AWS drive on DNAnexus:
```
dx api drive new '{"name": "VGP-test",
                   "cloud": "aws", 
                   "credentials": {"accessKeyId": [aws_access_key],  
                                   "secretAccessKey": [aws_secret_key]}}'
```

return:
```{
    "id": "drive-Kq7XjpkgZb9PzQYQGxBz8VvJ"
}```

To pull up user drives on DNAnexus:
```
dx api system findDrives '{"name": "VGP-test"}'
```
return:
```
{
    "results": [
        {
            "id": "drive-Kq7XjpkgZb9PzQYQGxBz8VvJ"
        }
    ]
}
```

To create symlink to existing AWS data object:
```
dx api file new '{"project": "project-xxx", 
                   "name": "filename",
                   "drive": "drive-xxx",
                   "symlinkPath": 
                   		{
                     	"container": "bucketname", 
                     	"object": "foldername/objectname" 
                   		}
                 }'
```
## Future Improvements
* Provide a mechanism for detecting files in AWS S3 which do not have DNAnexus links and creating them [ One option for this would be to add the DNAnexus file-id as metadata in AWS, and iterate through S3 files for objects that do not have this metadata field supplied. ]
* Provide mechanism for deleting DNAnexus links of AWS file is deleted
* Add enforcement or sanity checking for species ID naming convention
* Add mechanism for providing metadata to files during upload