package org.rhesusbase;

import java.awt.*;
import java.awt.List;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.*;

import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.Outlier;
import org.jfree.chart.renderer.OutlierList;
import org.jfree.chart.renderer.OutlierListCollection;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.CategoryItemRendererState;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.ui.RectangleEdge;

/**
 * Created with IntelliJ IDEA.
 * User: TheOne
 * Date: 15-11-27
 * Time: 下午11:17
 * To change this template use File | Settings | File Templates.
 */
public class boxplot extends BoxAndWhiskerRenderer {

    private String synonmousValue;
    private String nonsynonmousValue;
    private String cdsValue;
    private String utrValue;
    private String exonValue;
    private String intronValue;
    private String intergenicValue;
    private transient Paint medianPaint = Color.black;
    private transient Paint valuePaint = Color.red;
    private transient Paint boxPaint = Color.gray;

    public void setSynonmousValue(String value){
        this.synonmousValue = value;
    }

    public void setNonsynonmousValue (String value){
        this.nonsynonmousValue = value;
    }

    public void  setCdsValue (String value){
        this.cdsValue = value;
    }

    public void setUtrValue (String value){
        this.utrValue = value;
    }

    public void setExonValue (String value){
        this.exonValue = value;
    }

    public void setIntronValue (String value){
        this.intronValue = value;
    }

    public void setIntergenicValue (String value){
        this.intergenicValue = value;
    }

    public void setMedianPaint(Paint paint) {
        this.medianPaint = paint;
    }

    public void setValuePaint(Paint paint){
        this.valuePaint = paint;
    }

    public void setBoxPaint(Paint paint){
        this.boxPaint = paint;
    }

    public void drawHorizontalItem(Graphics2D g2,
                                   CategoryItemRendererState state, Rectangle2D dataArea,
                                   CategoryPlot plot, CategoryAxis domainAxis, ValueAxis rangeAxis,
                                   CategoryDataset dataset, int row, int column) {
        BoxAndWhiskerCategoryDataset bawDataset
                = (BoxAndWhiskerCategoryDataset) dataset;

        double categoryEnd = domainAxis.getCategoryEnd(column,
                getColumnCount(), dataArea, plot.getDomainAxisEdge());
        double categoryStart = domainAxis.getCategoryStart(column,
                getColumnCount(), dataArea, plot.getDomainAxisEdge());
        double categoryWidth = Math.abs(categoryEnd - categoryStart);

        double yy = categoryStart;
        int seriesCount = getRowCount();
        int categoryCount = getColumnCount();

        if (seriesCount > 1) {
            double seriesGap = dataArea.getHeight() * getItemMargin()
                    / (categoryCount * (seriesCount - 1));
            double usedWidth = (state.getBarWidth() * seriesCount)
                    + (seriesGap * (seriesCount - 1));
            // offset the start of the boxes if the total width used is smaller
            // than the category width
            double offset = (categoryWidth - usedWidth) / 2;
            yy = yy + offset + (row * (state.getBarWidth() + seriesGap));
        }
        else {
            // offset the start of the box if the box width is smaller than
            // the category width
            double offset = (categoryWidth - state.getBarWidth()) / 2;
            yy = yy + offset;
        }

        g2.setPaint(getItemPaint(row, column));
        Stroke s = getItemStroke(row, column);
        g2.setStroke(s);

        RectangleEdge location = plot.getRangeAxisEdge();

        Number xQ1 = bawDataset.getQ1Value(row, column);
        Number xQ3 = bawDataset.getQ3Value(row, column);
        Number xMax = bawDataset.getMaxRegularValue(row, column);
        Number xMin = bawDataset.getMinRegularValue(row, column);

        Shape box = null;
        if (xQ1 != null && xQ3 != null && xMax != null && xMin != null) {

            double xxQ1 = rangeAxis.valueToJava2D(xQ1.doubleValue(), dataArea,
                    location);
            double xxQ3 = rangeAxis.valueToJava2D(xQ3.doubleValue(), dataArea,
                    location);
            double xxMax = rangeAxis.valueToJava2D(xMax.doubleValue(), dataArea,
                    location);
            double xxMin = rangeAxis.valueToJava2D(xMin.doubleValue(), dataArea,
                    location);
            double yymid = yy + state.getBarWidth() / 2.0;
            double halfW = (state.getBarWidth() / 2.0) * this.getWhiskerWidth();

            // draw the box...
            box = new Rectangle2D.Double(Math.min(xxQ1, xxQ3), yy,
                    Math.abs(xxQ1 - xxQ3), state.getBarWidth());
            if (this.getFillBox()) {
                g2.setPaint(this.boxPaint);
                g2.fill(box);
            }

            Paint outlinePaint = getItemOutlinePaint(row, column);
            //if (this.useOutlinePaintForWhiskers) {
            //    g2.setPaint(outlinePaint);
            //}
            // draw the upper shadow...
            g2.draw(new Line2D.Double(xxMax, yymid, xxQ3, yymid));
            g2.draw(new Line2D.Double(xxMax, yymid - halfW, xxMax,
                    yymid + halfW));

            // draw the lower shadow...
            g2.draw(new Line2D.Double(xxMin, yymid, xxQ1, yymid));
            g2.draw(new Line2D.Double(xxMin, yymid - halfW, xxMin,
                    yy + halfW));

            g2.setStroke(getItemOutlineStroke(row, column));
            g2.setPaint(outlinePaint);
            g2.draw(box);
        }
        g2.setPaint(this.getArtifactPaint());

        g2.setPaint(this.valuePaint);
        double aRadius;                 // average radius
        //Number xMean = bawDataset.getMeanValue(row, column);
        //Number xMean = this.synonmousValue;
        Number xValue;
        if(column == 0 && this.synonmousValue != null){
            xValue = Double.parseDouble(this.synonmousValue);
        }else if (column == 1 && this.nonsynonmousValue != null){
            xValue = Double.parseDouble(this.nonsynonmousValue);
        }else if (column == 2 && this.cdsValue != null){
            xValue = Double.parseDouble(this.cdsValue);
        }else if (column == 3 && this.utrValue != null){
            xValue = Double.parseDouble(this.utrValue);
        }else if (column == 4 && this.exonValue != null){
            xValue = Double.parseDouble(this.exonValue);
        }else if (column == 5 && this.intronValue != null){
            xValue = Double.parseDouble(this.intronValue);
        }else if (column == 6 && this.intergenicValue != null){
            xValue = Double.parseDouble(this.intergenicValue);
        }else{
            xValue = null;
        }
        if (xValue != null) {
            double xxValue = rangeAxis.valueToJava2D(xValue.doubleValue(),
                    dataArea, location);
            aRadius = state.getBarWidth() / 8;
            // here we check that the average marker will in fact be
            // visible before drawing it...
            if ((xxValue > (dataArea.getMinX() - aRadius))
                    && (xxValue < (dataArea.getMaxX() + aRadius))) {
                Ellipse2D.Double avgEllipse = new Ellipse2D.Double(xxValue
                        - aRadius * 3.5, yy + aRadius * 0.5, aRadius * 1, aRadius * 1);
                g2.fill(avgEllipse);
                g2.draw(avgEllipse);
            }
        }

        g2.setPaint(this.medianPaint);
        if (this.isMedianVisible()) {
            Number xMedian = bawDataset.getMedianValue(row, column);
            if (xMedian != null) {
                double xxMedian = rangeAxis.valueToJava2D(xMedian.doubleValue(),
                        dataArea, location);
                g2.draw(new Line2D.Double(xxMedian, yy, xxMedian,
                        yy + state.getBarWidth()));
            }
        }

        if (state.getInfo() != null && box != null) {
            EntityCollection entities = state.getEntityCollection();
            if (entities != null) {
                addItemEntity(entities, dataset, row, column, box);
            }
        }
    }
    public void drawVerticalItem(Graphics2D g2, CategoryItemRendererState state,
                                 Rectangle2D dataArea, CategoryPlot plot, CategoryAxis domainAxis,
                                 ValueAxis rangeAxis, CategoryDataset dataset, int row, int column) {
        BoxAndWhiskerCategoryDataset bawDataset
                = (BoxAndWhiskerCategoryDataset) dataset;

        double categoryEnd = domainAxis.getCategoryEnd(column,
                getColumnCount(), dataArea, plot.getDomainAxisEdge());
        double categoryStart = domainAxis.getCategoryStart(column,
                getColumnCount(), dataArea, plot.getDomainAxisEdge());
        double categoryWidth = categoryEnd - categoryStart;

        double xx = categoryStart;
        int seriesCount = getRowCount();
        int categoryCount = getColumnCount();

        if (seriesCount > 1) {
            double seriesGap = dataArea.getWidth() * getItemMargin()
                    / (categoryCount * (seriesCount - 1));
            double usedWidth = (state.getBarWidth() * seriesCount)
                    + (seriesGap * (seriesCount - 1));
            // offset the start of the boxes if the total width used is smaller
            // than the category width
            double offset = (categoryWidth - usedWidth) / 2;
            xx = xx + offset + (row * (state.getBarWidth() + seriesGap));
        }
        else {
            // offset the start of the box if the box width is smaller than the
            // category width
            double offset = (categoryWidth - state.getBarWidth()) / 2;
            xx = xx + offset;
        }

        double yyAverage;
        double yyOutlier;

        Paint itemPaint = getItemPaint(row, column);
        g2.setPaint(itemPaint);
        Stroke s = getItemStroke(row, column);
        g2.setStroke(s);

        double aRadius = 0;                 // average radius

        RectangleEdge location = plot.getRangeAxisEdge();

        Number yQ1 = bawDataset.getQ1Value(row, column);
        Number yQ3 = bawDataset.getQ3Value(row, column);
        Number yMax = bawDataset.getMaxRegularValue(row, column);
        Number yMin = bawDataset.getMinRegularValue(row, column);
        Shape box = null;
        if (yQ1 != null && yQ3 != null && yMax != null && yMin != null) {

            double yyQ1 = rangeAxis.valueToJava2D(yQ1.doubleValue(), dataArea,
                    location);
            double yyQ3 = rangeAxis.valueToJava2D(yQ3.doubleValue(), dataArea,
                    location);
            double yyMax = rangeAxis.valueToJava2D(yMax.doubleValue(),
                    dataArea, location);
            double yyMin = rangeAxis.valueToJava2D(yMin.doubleValue(),
                    dataArea, location);
            double xxmid = xx + state.getBarWidth() / 2.0;
            double halfW = (state.getBarWidth() / 2.0) * this.getWhiskerWidth();

            // draw the body...
            box = new Rectangle2D.Double(xx, Math.min(yyQ1, yyQ3),
                    state.getBarWidth(), Math.abs(yyQ1 - yyQ3));
            Paint outlinePaint = getItemOutlinePaint(row, column);
            //if (this.useOutlinePaintForWhiskers) {
            //    g2.setPaint(outlinePaint);
            //}
            // draw the upper shadow...
            g2.draw(new Line2D.Double(xxmid, yyMax, xxmid, yyQ3));
            g2.draw(new Line2D.Double(xxmid - halfW, yyMax, xxmid + halfW, yyMax));

            // draw the lower shadow...
            g2.draw(new Line2D.Double(xxmid, yyMin, xxmid, yyQ1));
            g2.draw(new Line2D.Double(xxmid - halfW, yyMin, xxmid + halfW, yyMin));

            g2.setStroke(getItemOutlineStroke(row, column));
            g2.setPaint(outlinePaint);
            g2.draw(box);
        }

        g2.setPaint(this.getArtifactPaint());
        if (state.getInfo() != null && box != null) {
            EntityCollection entities = state.getEntityCollection();
            if (entities != null) {
                addItemEntity(entities, dataset, row, column, box);
            }
        }

        if (this.getFillBox()) {
            g2.setPaint(this.boxPaint);
            g2.fill(box);
        }

        g2.setPaint(this.medianPaint);
        if (this.isMedianVisible()) {
            Number yMedian = bawDataset.getMedianValue(row, column);
            if (yMedian != null) {
                double yyMedian = rangeAxis.valueToJava2D(
                        yMedian.doubleValue(), dataArea, location);
                g2.draw(new Line2D.Double(xx, yyMedian,
                        xx + state.getBarWidth(), yyMedian));
            }
        }

        //Number yMean = bawDataset.getMeanValue(row, column);
        //Number yMean = this.synonmousValue;
        g2.setPaint(this.valuePaint);
        Number yValue;
        if(column == 0 && this.synonmousValue != null){
            yValue = Double.parseDouble(this.synonmousValue);
        }else if (column == 1 && this.nonsynonmousValue != null){
            yValue = Double.parseDouble(this.nonsynonmousValue);
        }else if (column == 2 && this.cdsValue != null){
            yValue = Double.parseDouble(this.cdsValue);
        }else if (column == 3 && this.utrValue != null){
            yValue = Double.parseDouble(this.utrValue);
        }else if (column == 4 && this.exonValue != null){
            yValue = Double.parseDouble(this.exonValue);
        }else if (column == 5 && this.intronValue != null){
            yValue = Double.parseDouble(this.intronValue);
        }else if (column == 6 && this.intergenicValue != null){
            yValue = Double.parseDouble(this.intergenicValue);
        }else{
            yValue = null;
        }
        if (yValue != null) {
            yyAverage = rangeAxis.valueToJava2D(yValue.doubleValue(),
                    dataArea, location);
            aRadius = state.getBarWidth() / 8;
            // here we check that the average marker will in fact be
            // visible before drawing it...
            if ((yyAverage > (dataArea.getMinY() - aRadius))
                    && (yyAverage < (dataArea.getMaxY() + aRadius))) {
                Ellipse2D.Double avgEllipse = new Ellipse2D.Double(
                        xx + aRadius * 3.5, yyAverage - aRadius * 0.5, aRadius * 1,
                        aRadius * 1);
                g2.fill(avgEllipse);
                g2.draw(avgEllipse);
            }
        }
    }
}
