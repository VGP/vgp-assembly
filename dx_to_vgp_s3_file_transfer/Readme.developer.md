# Platform to AWS s3 Developer Readme

## Notes for developers
AWS recommends not uploading files with initial common paths i.e. 

>**p**roject\itemcontain\file1  
>**p**roject\anotheritem\file2  
>**p**roject\thirditem\file3  


The common '**p**' start to all the paths will result in future [performance issues with S3 look-ups](https://aws.amazon.com/blogs/aws/amazon-s3-performance-tips-tricks-seattle-hiring-event/) when querying and making S3 request.

## App timeout

This app has a hard timeout set of 48hrs. If a transfer will take longer than 48hrs, feel free to alter the `runSpec.timeoutPolicy` field in the *dxapp.json*.


## Running this app with additional computational resources

This app has the following entry points:

* main
* create_upload_report
* s3_upload

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See <a
href="https://wiki.dnanexus.com/API-Specification-v1.0.0/IO-and-Run-Specifications#Run-Specification">Run
Specification</a> in the API documentation for more information about the
available instance types.
