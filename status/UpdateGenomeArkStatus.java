import java.util.HashMap;

import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.S3ObjectSummary;

public class UpdateGenomeArkStatus {
	
	private static final int SPECIES_NAME = 1;
	private static final int SPECIES_ID = 2;
	private static final int DATA_TYPE = 3;
	private static final int SEQ_PLATFORM = 4;
	private static final int FILE = 5;

	public static void main(String[] args) {
		
		if (args.length > 0 && args[0].contains("help")) {
			System.out.println("Usage: java -jar UpdateGenomeArkStatus.jar [-md]");
			System.exit(0);
		}
		
		boolean isMDstyle = false;
		if (args.length > 0 && args[0].contains("md")) {
			isMDstyle = true;
		}
		
        AmazonS3 s3 = AmazonS3ClientBuilder.standard()
            .withRegion("us-east-1")
            .build();

        String bucketName = "genomeark";

        try {
        	/*
             * List all objects under s3://genomeark/species
             * and add the datatype.
             */
            
            ListObjectsRequest listObjectsRequest = new ListObjectsRequest()
                    .withBucketName(bucketName)
                    .withPrefix("species");
            ObjectListing objectListing = s3.listObjects(listObjectsRequest);
            Double countKeys = 0d;
            
            String[] tokens;
            String speciesId;
            String seqPlatform;
            String file;
            HashMap<String, VGPData> dataMap = new HashMap<String, VGPData>();
            VGPData data;
            while (true) {
            	for (S3ObjectSummary objectSummary : objectListing.getObjectSummaries()) {
            		//System.err.println(objectSummary.getKey());
            		countKeys++;
            		tokens = objectSummary.getKey().split("/");
            		if (tokens.length < 6) continue;
            		if (tokens[DATA_TYPE].equals("genomic_data")) {
            			speciesId = tokens[SPECIES_ID];
            			if (dataMap.containsKey(speciesId)) {
            				data = dataMap.get(speciesId);
            			} else {
            				data = new VGPData(speciesId, tokens[SPECIES_NAME]);
            				dataMap.put(speciesId, data);
            			}
            			seqPlatform = tokens[SEQ_PLATFORM].toLowerCase();
            			file = tokens[FILE].toLowerCase();
            			if (seqPlatform.equals("pacbio")) {
            				data.addPacbio(file);
            			} else if (seqPlatform.equals("10x")) {
            				data.add10xR1(file);
            			} else if (seqPlatform.equals("bionano")) {
            				data.addBionano(file);
            			} else {
            				data.addHiC(seqPlatform);
            			}
            		}
            	}
            	if (objectListing.isTruncated()) {
            		objectListing = s3.listNextBatchOfObjects(objectListing);
            	} else {
            		break;
            	}
            }
            System.err.println("Total num. objects under s3://genomeark/species: " + String.format("%,.0f", countKeys));
            
            VGPData.printHeader(isMDstyle);
            for (String id : dataMap.keySet()) {
            	dataMap.get(id).printVGPData(isMDstyle);
            }
        } catch (AmazonServiceException ase) {
            System.err.println("Caught an AmazonServiceException, which means your request made it "
                    + "to Amazon S3, but was rejected with an error response for some reason.");
            System.err.println("Error Message:    " + ase.getMessage());
            System.err.println("HTTP Status Code: " + ase.getStatusCode());
            System.err.println("AWS Error Code:   " + ase.getErrorCode());
            System.err.println("Error Type:       " + ase.getErrorType());
            System.err.println("Request ID:       " + ase.getRequestId());
        } catch (AmazonClientException ace) {
            System.err.println("Caught an AmazonClientException, which means the client encountered "
                    + "a serious internal problem while trying to communicate with S3, "
                    + "such as not being able to access the network.");
            System.err.println("Error Message: " + ace.getMessage());
        }
	}

}
