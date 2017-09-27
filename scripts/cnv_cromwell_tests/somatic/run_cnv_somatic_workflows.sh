#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/cnv_common_tasks.wdl

WORKING_DIR=/home/travis/build/broadinstitute

pushd .
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk
# IMPORTANT: This code is duplicated in the M2 WDL test.
if [ ${TRAVIS_PULL_REQUEST} != false ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${TRAVIS_PULL_REQUEST};
else
  HASH_TO_USE=${TRAVIS_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
fi
echo "Docker build done =========="

popd

echo "Inserting docker image into json ========"
CNV_CROMWELL_TEST_DIR="${WORKING_DIR}/gatk/scripts/cnv_cromwell_tests/somatic/"
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wes_workflow.json >cnv_somatic_panel_wes_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wgs_workflow.json >cnv_somatic_panel_wgs_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wes_do-gc_workflow.json >cnv_somatic_panel_wes_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wgs_do-gc_workflow.json >cnv_somatic_panel_wgs_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_bam_wes_workflow.json >cnv_somatic_bam_wes_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_bam_wgs_workflow.json >cnv_somatic_bam_wgs_workflow_mod.json

echo "Running ========"
CROMWELL_JAR="cromwell-0.28.jar"

# Panel WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wes_workflow_mod.json
# Panel WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wgs_workflow_mod.json
# Panel WES w/ explicit GC correction
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wes_do-gc_workflow_mod.json
# Panel WGS w/ explicit GC correction
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wgs_do-gc_workflow_mod.json

# BAM WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_bam_workflow.wdl cnv_somatic_bam_wes_workflow_mod.json
# BAM WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_bam_workflow.wdl cnv_somatic_bam_wgs_workflow_mod.json