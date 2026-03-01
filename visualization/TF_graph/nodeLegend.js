//import * as d3 from 'd3';

// based on https://vizhub.com/curran/a446f43c024a49608f7ae418cde946a2?file=colorLegend.js
export const nodeLegend = (
  selection,
  {
    colorScale,
    colorLabel,
    sizeScale,
    sizeLabel,
    legendLabelFontSize,
    legendX,
    legendY,
    legendFontSize,
    legendPadding,
    colorPadding,
    colorHeight,
    colorWidth,
    colorTickNumber,
    colorNumber,
    legendSpacing,
    sizeNumber,
    sizeSpacing,
    sizePadding,
  },
) => {
  const colorGroup = selection
    .selectAll('g.color-legend')
    .data([null])
    .join('g')
    .attr('class', 'color-legend')
    .attr('transform', `translate(${legendX},${legendY})`);

  colorGroup
    .selectAll('text.label')
    .data([null])
    .join('text')
    .attr('class', 'label')
    .attr('font-size', legendLabelFontSize)
    .text(colorLabel);

  const colorTicks = colorScale.ticks(colorTickNumber);

  colorGroup
    .selectAll('text.ticks')
    .data(colorTicks)
    .join('text')
    .attr('class', 'ticks')
    .text((d) => d)
    .attr('x', colorWidth + 5)
    .attr(
      'y',
      (d, i) =>
        legendPadding +
        (i * colorHeight) / colorTicks.length,
    )
    .attr('dy', '0.32em')
    .attr('font-size', legendFontSize);
  //.attr('dominant-baseline', 'central')
  //.attr('text-alignment', 'begin');

  const legendColorVals = colorScale.ticks(colorNumber);

  colorGroup
    .selectAll('rect.heat')
    .data(legendColorVals)
    .join('rect')
    .attr('class', 'heat')
    .attr('x', 0)
    .attr(
      'y',
      (d, i) =>
        colorPadding +
        ((i - 0.5) * colorHeight) / legendColorVals.length,
    )
    .attr('height', colorHeight / legendColorVals.length)
    .attr('width', colorWidth)
    .attr('fill', (d) => colorScale(d));

  const sizeGroup = selection
    .selectAll('g.size-legend')
    .data([null])
    .join('g')
    .attr('class', 'size-legend')
    .attr(
      'transform',
      `translate(${legendX},${legendY + colorHeight + legendSpacing})`,
    );

  sizeGroup
    .selectAll('text.size-label')
    .data([null])
    .join('text')
    .attr('class', 'size-label')
    .attr('font-size', legendLabelFontSize)
    .text(sizeLabel);

  const sizeTicks = sizeScale.ticks(sizeNumber);

  sizeGroup
    .selectAll('g.size-ticks')
    .data(sizeTicks)
    .join((enter) =>
      enter
        .append('g')
        .attr('class', 'tick')
        .call((selection) => {
          selection.append('text');
          selection.append('circle');
        }),
    )
    .attr('class', 'size-ticks')
    .attr(
      'transform',
      (d, i) =>
        `translate(15, ${legendPadding + i * sizeSpacing})`,
    )
    .call((selection) => {
      selection
        .select('circle')
        .attr('cx', 0)
        .attr('cy', 0)
        .attr('r', (d) => sizeScale(d));

      selection
        .select('text')
        .text((d) => d)
        .attr('x', sizePadding)
        .attr('y', 0)
        .attr('dominant-baseline', 'middle')
        .attr('font-size', legendFontSize);
    });
};
