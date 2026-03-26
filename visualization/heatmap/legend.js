import * as d3 from 'd3';

// based on https://vizhub.com/curran/a446f43c024a49608f7ae418cde946a2?file=colorLegend.js
export const legend = (
  svg,
  {
    colorScale,
    legendLabel,
    legendX,
    legendY,
    legendPadding,
    legendHeight,
    legendWidth,
    legendTickNumber,
    legendColorNumber,
  },
) => {
  const legendGroup = svg
    .selectAll('g.legend')
    .data([null])
    .join('g')
    .attr('class', 'legend')
    .attr('transform', `translate(${legendX},${legendY})`);

  legendGroup
    .selectAll('text.legend-label')
    .data([null])
    .join('text')
    .attr('class', 'legend-label')
    .text(legendLabel);

  const legendTicks = colorScale.ticks(legendTickNumber);

  legendGroup
    .selectAll('text.legend-ticks')
    .data(legendTicks)
    .join('text')
    .attr('class', 'legend-ticks')
    .text((d) => d)
    .attr('x', legendWidth + 5)
    .attr(
      'y',
      (d, i) =>
        legendPadding +
        (i * legendHeight) / legendTicks.length,
    )
        .attr('dy', '0.32em')
    //.attr('dominant-baseline', 'central')
    //.attr('text-alignment', 'begin');

  const legendColorVals = colorScale.ticks(
    legendColorNumber,
  );

  legendGroup
    .selectAll('rect.legend-heat')
    .data(legendColorVals)
    .join('rect')
    .attr('class', 'legend-heat')
    .attr('x', 0)
    .attr(
      'y',
      (d, i) =>
        legendPadding +
        ((i - 0.5) * legendHeight) / legendColorVals.length,
    )
    .attr('height', legendHeight / legendColorVals.length)
    .attr('width', legendWidth)
    .attr('fill', (d) => colorScale(d));
};
