# How to upload a file to the AWS GenomeArk

There are many ways to upload your data.

## Upload from a local cluster

### Obtain credentials
Contact the Genome10K VGP (vgp-assembly @ googlegroups.com) to obtain access credentials for uploading.
Of note, the VGP does not require any credentials for downloading. By downloading the data,
you agree and accept the [data use policy](https://genome10k.soe.ucsc.edu/about/data_use_policy#embargo).

### Set-up AWS CLI
Install AWS CLI (Command Line Interface) from [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html).
Bellow is a short example if you *don't* have root permission.

```
curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip”
unzip awscli-bundle.zip
./awscli-bundle/install -i $path_to_install
export PATH=$path_to_install/bin:$PATH
```

The next step is to configure your credentials obtained from the assembly-group.
```
aws configure
```
Type in your given ws_access_key_id and aws_secret_access_key_id when prompted.

Copy a file following the [data_structure]("https://github.com/VGP/vgp-assembly/blob/master/DNAnexus_and_AWS_data_structure.md").
```
aws s3 cp <file> s3://genomeark-upload/species/<species_name>/<species_id>/<data_type>/<file>
```
For example, uploading a pacbio subread.bam from the hummingbird will be
```
aws s3 mXXXX.subreads.bam s3://genomeark-upload/species/Calypte_anna/bCalAnn1/pacbio/mXXXX.subreads.bam
```

### Check your file after transfer is completed
Run [check_etag.sh](vgp-assembly/aws_upload/utils/check_etag.sh) and see if it matches the eTag on the uploaded file.
The eTag will be the md5 (or md5sum) for files <5 Gb, and a combined hash of multi-part files when larger than 5Gb.
```
./check_etag.sh mXXXX.subreads.bam
```
Here is a MacOS version of check_etag.sh: 

There are ways to change the eTag, but please use the default behavior of the aws cli and not change the eTag.

### Contact us back
This is *very* important. Let us know when your uploading is completed.
After a short check for the file structure, your data will be transferred to the GenomeArk.

## Transfer from DNAnexus directly to the GenomeArk
TBA
