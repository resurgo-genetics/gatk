package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.rsvd.HDF5RandomizedSVDReadCountPanelOfNormals.HDF5ChunkedDoubleMatrixHelper;
import org.broadinstitute.hellbender.tools.pon.PoNTestUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * TODO
 */
public final class HDF5RandomizedSVDReadCountPanelOfNormalsUnitTest extends BaseTest {
    private static final int CHUNK_DIVISOR = 256;

    @DataProvider(name = "testCreateLargeMatrixData")
    public Object[][] dataProvider() {
        return new Object[][] {
                new Object[] {100, 100000},
                new Object[] {CHUNK_DIVISOR + 16, (Integer.MAX_VALUE / Byte.SIZE) / CHUNK_DIVISOR},
                new Object[] {(Integer.MAX_VALUE / Byte.SIZE) / CHUNK_DIVISOR, CHUNK_DIVISOR + 16}
        };
    }

    @Test(dataProvider = "testCreateLargeMatrixData")
    public void testCreateLargeMatrix(final int numRows,
                                      final int numColumns) {
        final double mean = 0.;
        final double sigma = 1.;
        final String matrixPath = "/test/matrix";

        final RealMatrix largeMatrix = createMatrixOfGaussianValues(numRows, numColumns, mean, sigma);
        logger.debug("Large matrix created.");
        final File tempOutputHD5 = IOUtils.createTempFile("large-matrix-", ".hd5");
        try (final HDF5File hdf5File = new HDF5File(tempOutputHD5, HDF5File.OpenMode.CREATE)) {
            HDF5ChunkedDoubleMatrixHelper.writeMatrix(hdf5File, matrixPath, largeMatrix.getData(), CHUNK_DIVISOR);
        }
        logger.debug("Large matrix written.");

        try (final HDF5File hdf5FileForReading = new HDF5File(tempOutputHD5, HDF5File.OpenMode.READ_ONLY)) {
            final double[][] result = HDF5ChunkedDoubleMatrixHelper.readMatrix(hdf5FileForReading, matrixPath);
            final RealMatrix resultAsRealMatrix = new Array2DRowRealMatrix(result, false);
            Assert.assertTrue(resultAsRealMatrix.getRowDimension() == numRows);
            Assert.assertTrue(resultAsRealMatrix.getColumnDimension() == numColumns);
            PoNTestUtils.assertEqualsMatrix(resultAsRealMatrix, largeMatrix, false);
        }
    }

    private static RealMatrix createMatrixOfGaussianValues(final int numRows,
                                                           final int numColumns,
                                                           final double mean,
                                                           final double sigma) {
        final RealMatrix largeMatrix = new Array2DRowRealMatrix(numRows, numColumns);
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        randomDataGenerator.reSeed(337337337);
        largeMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return randomDataGenerator.nextGaussian(mean, sigma);
            }
        });
        return largeMatrix;
    }
}