import { legend } from './legend.js';

export const body = (
  selection,
  {
    bodyData,
    organDict,
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
  const outlineData = bodyData.filter(
    (d) => d.id === 'body',
  );

  const organData = bodyData.filter((d) => d.id !== 'body');

  selection
    .selectAll('path.outline')
    .data(outlineData)
    .join('path')
    .attr('class', 'outline')
    .attr('fill', (d) => d.fill)
    .attr('stroke', (d) => d.stroke)
    .attr('stroke-width', (d) =>
      d.strokeWidth ? d.strokeWidth : 1,
    )
    .attr('d', (d) => d.d)
    .attr('transform', (d) =>
      d.rotate
        ? `translate(${d.x}, ${d.y}) scale(${d.scale}) rotate(${d.rotate})`
        : `translate(${d.x}, ${d.y}) scale(${d.scale})`,
    );

  selection
    .selectAll('path.organ')
    .data(organData)
    .join('path')
    .attr('class', 'organ')
    .attr('fill', (d) => d.fill)
    .attr('stroke', (d) => d.stroke)
    .attr('stroke-width', (d) =>
      d.strokeWidth ? d.strokeWidth : 1,
    )
    .attr('d', (d) => d.d)
    .attr('transform', (d) =>
      d.rotate
        ? `translate(${d.x}, ${d.y}) scale(${d.scale}) rotate(${d.rotate})`
        : `translate(${d.x}, ${d.y}) scale(${d.scale})`,
    )
    .attr('opacity', (d) =>
      clickedOrgan1 || hoveredOrgan
        ? d.id === organDict[clickedOrgan1] ||
          d.id === organDict[clickedOrgan2] ||
          d.id === organDict[hoveredOrgan]
          ? 1
          : 0.2
        : 1,
    );

  legend(selection, {
    dataDict: organDict,
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
  });
};
