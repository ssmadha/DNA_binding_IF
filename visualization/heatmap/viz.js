
export const viz = (svg, { state, setState }) => {
  const {
    data,
    xAnnotation,
    yAnnotation,
    yAnnoLabels,
    xLabels,
    xRotation,
    xAxisLabelText,
    xAxisLabelOffset,
    yLabels,
    yTickFormat = '',
    yTickSize,
    yAxisLabelText,
    yAxisLabelOffset,
    heatValue,
    marginLeft,
    marginTop,
    marginRight,
    marginBottom,
    width,
    height,
    innerRectFill,
    padding,
    xAnnoPadding,
    xAnnoHeight,
    yAnnoPadding,
    yAnnoWidth,
    yFilter,
    legendLabel,
    legendX,
    legendY,
    legendPadding,
    legendHeight,
    legendWidth,
    legendTickNumber,
    legendColorNumber,
  } = state;

  const filtData =
    yFilter !== undefined
      ? data.filter(
          (d) => yAnnotation[yAnnoLabels(d)] === yFilter,
        )
      : data;

  const xScale = d3
    .scaleBand(filtData.map(xLabels), [
      marginLeft,
      width - marginRight,
    ])
    .padding(padding);

  const yScale = d3.scaleBand(filtData.map(yLabels), [
    height - marginBottom,
    marginTop,
  ]);

  const heatScale = d3.scaleSequential(
    [0, d3.max(filtData, heatValue)],
    d3.interpolateInferno,
  );

  const xColor = d3.scaleOrdinal(
    Object.values(xAnnotation),
    d3.schemePaired,
  );

  const yColor = d3.scaleOrdinal(
    Object.values(yAnnotation),
    d3.schemePaired,
  );

  const innerWidth = width - marginLeft - marginRight;
  const innerHeight = height - marginTop - marginBottom;

  svg
    .selectAll('rect.inner-rectangle')
    .data([null])
    .join('rect')
    .attr('class', 'inner-rectangle')
    .attr('x', marginLeft)
    .attr('y', marginTop)
    .attr('width', innerWidth)
    .attr('height', innerHeight)
    .attr('fill', innerRectFill);

  axes(svg, {
    height,
    width,
    marginBottom,
    marginLeft,
    marginRight,
    marginTop,
    xScale,
    xRotation,
    yScale,
    yTickFormat,
    yTickSize,
    xAxisLabelText,
    xAxisLabelOffset,
    yAxisLabelText,
    yAxisLabelOffset,
  });

  const hmToolTip = d3
    .select('body')
    .append('div')
    .attr('class', 'tooltip')
    .style('opacity', 0);

  // Heatmap
  svg
    .selectAll('rect.mark')
    .data(filtData)
    .join('rect')
    .attr('class', 'mark')
    .attr('x', (d) => xScale(xLabels(d)))
    .attr('y', (d) => yScale(yLabels(d)))
    .attr('width', xScale.bandwidth())
    .attr('height', yScale.bandwidth())
    .attr('fill', (d) => heatScale(heatValue(d)))
    .on('mouseover', function (event, d) {
      const [pointerX, pointerY] = d3.pointer(event);
      hmToolTip
        .transition()
        .duration(200)
        .style('opacity', 0.9);
      hmToolTip
        .html(
          `Gene: ${yLabels(d)} <br>
        Tissue: ${xLabels(d)}`,
        )
        .style('position', 'fixed')
        .style('left', event.pageX + 10 + 'px')
        .style('top', event.pageY + 10 + 'px')
        .style('background', '#bbbbbb');
    })
    .on('mouseout', function () {
      hmToolTip
        .transition()
        .duration(500)
        .style('opacity', 0);
      hmToolTip.html('');
    });

  //remove duplicates from the xLabels column of filtData
  const xAnnotationDedup = Array.from(
    new Set(filtData.map(xLabels)),
  );

  //X annotation
  svg
    .selectAll('rect.xAnnotation')
    .data(xAnnotationDedup)
    .join('rect')
    .attr('class', 'xAnnotation')
    .attr('x', (d) => xScale(d))
    .attr('y', marginTop - xAnnoPadding - xAnnoHeight)
    .attr('width', xScale.bandwidth())
    .attr('height', xAnnoHeight)
    .attr('fill', (d) => xColor(xAnnotation[d]));

  //remove duplicates from the yLabels and yAnnoLabels columns of filtData
  const yAnnotationDedup = Array.from(
    new Set(
      filtData.map((d) => [yLabels(d), yAnnoLabels(d)]),
    ),
  );

  //Y annotation
  svg
    .selectAll('rect.yAnnotation')
    .data(yAnnotationDedup)
    .join('rect')
    .attr('class', 'yAnnotation')
    .attr('x', width - marginRight + yAnnoPadding)
    .attr('y', (d) => yScale(d[0]))
    .attr('width', yAnnoWidth)
    .attr('height', yScale.bandwidth())
    .attr('fill', (d) => yColor(yAnnotation[d[1]]))
    .style('cursor', 'pointer')
    .on('click', (event, d) => {
      setState((state) => ({
        ...state,
        yFilter:
          yFilter === undefined
            ? yAnnotation[d[1]]
            : undefined,
        yTickFormat: yFilter === undefined ? null : '',
      }));
    });

  legend(svg, {
    colorScale: heatScale,
    legendLabel,
    legendX,
    legendY,
    legendPadding,
    legendHeight,
    legendWidth,
    legendTickNumber,
    legendColorNumber,
  });
};
