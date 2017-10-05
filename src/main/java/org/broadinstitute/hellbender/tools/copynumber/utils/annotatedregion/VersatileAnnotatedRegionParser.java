package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.supercsv.io.CsvMapReader;
import org.supercsv.io.ICsvMapReader;
import org.supercsv.prefs.CsvPreference;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

// TODO: Catch exceptions and handle with GATK errors.

/**
 * Generally a good idea to create a new instance for each file being read.
 */
public class VersatileAnnotatedRegionParser {
    protected final static Set<String> CONTIG_HEADERS = Sets.newHashSet("Chromosome", "contig", "chr", "CONTIG");
    protected final static Set<String> START_HEADERS = Sets.newHashSet("Start_position", "Start_Position", "start", "START");
    protected final static Set<String> END_HEADERS = Sets.newHashSet("End_position", "End_Position", "end", "stop", "END");

    protected HashCommentMatcher hashCommentMatcher = new HashCommentMatcher();

    /**
     *  TODO: Docs
     * @param tsvRegionFile -- File containing tsv of genomic regions and annotations per line.  E.g. a segment file.
     * @param headersOfInterest -- should not include any headers that are used to define the region (e.g. contig, start, end)
     * @return
     * @throws IOException
     */
    public List<SimpleAnnotatedGenomicRegion> readAnnotatedRegions(final File tsvRegionFile, final Set<String> headersOfInterest) throws IOException {

        // TODO: Get rid of this line and make attribute final
        hashCommentMatcher = new HashCommentMatcher();

        // Tab-delimited text file with comments on lines starting with "#"
        final CsvPreference tsvPreference = new CsvPreference.Builder('"', '\t', "\n")
                .skipComments(hashCommentMatcher).build();

        ICsvMapReader mapReader = null;
        final List<SimpleAnnotatedGenomicRegion> result = new ArrayList<>();
        try {
            mapReader = new CsvMapReader(new FileReader(tsvRegionFile), tsvPreference);

            // the header columns are used as the keys to the Map
            final String[] headers = mapReader.getHeader(true);

            // Create mapping to get the contig, start, and end.
            final Set<String> headerSet = new HashSet<>(Arrays.asList(headers));
            final Set<String> contigHeaders = Sets.intersection(headerSet, CONTIG_HEADERS);
            final Set<String> startHeaders = Sets.intersection(headerSet, START_HEADERS);
            final Set<String> endHeaders = Sets.intersection(headerSet, END_HEADERS);

            validateInputHeader(contigHeaders, startHeaders, endHeaders);

            final Set<String> otherHeaders = Sets.intersection(headerSet, headersOfInterest);
            final String contigHeader = contigHeaders.iterator().next();
            final String startHeader = startHeaders.iterator().next();
            final String endHeader = endHeaders.iterator().next();

            Map<String, String> segmentInfoMap;
            while ((segmentInfoMap = mapReader.read(headers)) != null) {
                final Map<String, String> annotations = segmentInfoMap.entrySet().stream()
                        .filter(e -> otherHeaders.contains(e.getKey()))
                        .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue()));

                result.add(new SimpleAnnotatedGenomicRegion(
                        new SimpleInterval(segmentInfoMap.get(contigHeader),
                                Integer.parseInt(segmentInfoMap.get(startHeader)),
                                Integer.parseInt(segmentInfoMap.get(endHeader))),
                        annotations));
            }

        } finally {
            if (mapReader != null) {
                mapReader.close();
            }
        }

        return result;
    }

    private void validateInputHeader(Set<String> contigHeaders, Set<String> startHeaders, Set<String> endHeaders) {
        if (contigHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for contig.  Must be one.  Found: " + contigHeaders.toString());
        }
        if (startHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for start.  Must be one.  Found: " + startHeaders.toString());
        }
        if (endHeaders.size() != 1) {
            throw new UserException.BadInput("Wrong number of headers that could be used for end.  Must be one.  Found: " + endHeaders.toString());
        }
    }

    // TODO: Make unmodifiable?
    public List<String> getComments() {
        return hashCommentMatcher.getComments();
    }

    public void writeAnnotatedRegionsAsTsv(final List<SimpleAnnotatedGenomicRegion> regions, final File outputFile, final List<String> comments) {

    }
}
