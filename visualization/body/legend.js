//import * as d3 from 'd3';

// based on https://vizhub.com/curran/a446f43c024a49608f7ae418cde946a2?file=colorLegend.js
export const legend = (
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
    clickedOrgan1,
    clickedOrgan2,
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
      clickedOrgan1 || hoveredOrgan
        ? d === clickedOrgan1 ||
          d === clickedOrgan2 ||
          d === hoveredOrgan
          ? 1
          : 0.2
        : 1,
    )
    .on('click', (event, d) => {
      clickedOrgan1
        ? clickedOrgan1 === d
          ? setClickedOrgan(null)
          : clickedOrgan2
            ? clickedOrgan2 === d
              ? setClickedOrgan(clickedOrgan1, null)
              : setClickedOrgan(clickedOrgan1, d)
            : setClickedOrgan(clickedOrgan1, d)
        : setClickedOrgan(d);
    })
    .on('mouseover', (event, d) => setHoveredOrgan(d))
    .on('mouseout', () => setHoveredOrgan(null))
    .style('cursor', 'pointer');
};
