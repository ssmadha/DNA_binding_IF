//import * as d3 from 'd3';

// based on https://vizhub.com/curran/a446f43c024a49608f7ae418cde946a2?file=colorLegend.js
export const bodyLegend = (
  svg,
  {
    dataDict,
    legendLabel,
    legendLabelFontSize,
    legendX,
    legendY,
    legendPadding,
    legendFontSize,
    tickSpacing,
    tickPadding,
    clickedOrgan,
    setClickedOrgan,
    hoveredOrgan,
    setHoveredOrgan,
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
    .text(legendLabel)
    .attr('font-size', legendLabelFontSize);

  legendGroup
    .selectAll('g.tick')
    .data(Object.keys(dataDict))
    .join((enter) =>
      enter
        .append('g')
        .attr('class', 'tick')
        .call((selection) => {
          selection.append('text');
        }),
    )
    .attr(
      'transform',
      (d, i) =>
        `translate(0,${legendPadding + i * tickSpacing})`,
    )
    .call((selection) => {
      selection
        .select('text')
        .text((d) => d)
        .attr('dy', '0.32em')
        .attr('x', tickPadding);
    })
    .attr('font-size', legendFontSize)
    .attr('opacity', (d) =>
      clickedOrgan || hoveredOrgan
        ? d === clickedOrgan || d === hoveredOrgan
          ? 1
          : 0.2
        : 1,
    )
    .on('mouseenter', (event, d) => setHoveredOrgan(d))
    .on('mouseleave', () => setHoveredOrgan(null))
    .on('click', (event, d) => {
      clickedOrgan === d
        ? setClickedOrgan(null)
        : setClickedOrgan(d);
    })
    .style('cursor', 'pointer');
};
