#!/usr/bin/env bash

# Copy genome data to HDFS on a GCS cluster.

${GATK_HOME:-../..}/gatk-launch ParallelCopyGCSDirectoryIntoHDFSSpark \
    --inputGCSPath gs://broad-spark-eval-test-data/genome/ \
    --outputHDFSDirectory hdfs://${GCS_CLUSTER}-m:8020/user/$USER/q4_spark_eval \
    -apiKey $API_KEY \
    -- \
    --sparkRunner GCS \
    --cluster $GCS_CLUSTER