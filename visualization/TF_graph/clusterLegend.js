//import * as d3 from 'd3';

// based on https://vizhub.com/curran/a446f43c024a49608f7ae418cde946a2?file=colorLegend.js
export const clusterLegend = (
  selection,
  {
    colorScale,
    legendLabel,
    legendLabelY,
    legendPadding,
    tickSpacing,
    legendCircleRadius,
    tickPadding,
    hoveredCluster,
    setHoveredCluster,
    clickedCluster,
    setClickedCluster,
  },
) => {
  const legendGroup = selection
    .selectAll('g.legend')
    .data([null])
    .join('g')
    .attr('class', 'legend');

  legendGroup
    .selectAll('text.legend-label')
    .data([null])
    .join('text')
    .attr('class', 'legend-label')
    .attr('x', 0)
    .attr('y', legendLabelY)
    .text(legendLabel);

  legendGroup
    .selectAll('g.tick')
    .data(colorScale.domain())
    .join((enter) =>
      enter
        .append('g')
        .attr('class', 'tick')
        .call((selection) => {
          selection.append('text');
          selection.append('circle');
        }),
    )
    .attr(
      'transform',
      (d, i) =>
        `translate(${legendPadding + i * tickSpacing}, 15)`,
    )
    .call((selection) => {
      selection
        .select('text')
        .text((d) => d)
        // .attr('x', 0)
        // .attr('y', 0)
        .attr('text-anchor', 'middle');
      selection
        .select('circle')
        .attr('cx', 0)
        .attr('cy', tickPadding)
        .attr('r', legendCircleRadius)
        .attr('fill', 'white')
        .attr('stroke', colorScale)
        .attr('stroke-width', 2);
    })
    .attr('opacity', (d) =>
      (hoveredCluster !== null &&
        hoveredCluster !== undefined) ||
      (clickedCluster !== null &&
        clickedCluster !== undefined)
        ? d === hoveredCluster || d === clickedCluster
          ? 1
          : 0.2
        : 1,
    )
    .on('mouseover', (event, d) => {
      setHoveredCluster(d);
    })
    .on('mouseout', () => setHoveredCluster(null))
    .on('click', (event, d) => {
      if (d === clickedCluster) {
        setClickedCluster(null);
      } else {
        setClickedCluster(d);
      }
    })
    .style('cursor', 'pointer');
};
