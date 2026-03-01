//import * as d3 from 'd3';

export const graph = (
  selection,
  {
    nodes,
    links,
    packData,
    nodeStrokeWidth = 1,
    linkStrokeWidth = 1,
    textFontSize = 4,
    transform,
    zoomNameLimit,
    hoveredCluster,
    clickedCluster,
    setClickedNode,
    setState,
  },
) => {
  //drag code from https://vizhub.com/curran/a6c068082da148f3b637e7c5b4bfabcb
  const dragBehavior = d3
    .drag()
    .on('drag', (event, draggedDatum) => {
      setState((state) => ({
        ...state,
        nodes: state.nodes.map((d) =>
          d.id === draggedDatum.id
            ? { ...d, x: event.x, y: event.y }
            : d,
        ),
        links: state.links.map((d) =>
          d.source.id === draggedDatum.id
            ? {
                ...d,
                source: {
                  ...d.source,
                  x: event.x,
                  y: event.y,
                },
              }
            : d.target.id === draggedDatum.id
              ? {
                  ...d,
                  target: {
                    ...d.target,
                    x: event.x,
                    y: event.y,
                  },
                }
              : d,
        ),
        freezeSim: true,
      }));
    });

  let hmToolTip = d3.select('body .tooltip');

  if (hmToolTip.empty()) {
    hmToolTip = d3
      .select('body')
      .append('div')
      .attr('class', 'tooltip')
      .style('opacity', 0);
  }

  selection
    .selectAll('rect.background')
    .data([null])
    .join('rect')
    .attr('class', 'background')
    .attr('width', '100%')
    .attr('height', '100%')
    .attr('fill', '#FFFFFF');

  const networkGroup = selection
    .selectAll('g.networkGraph')
    .data([null])
    .join('g')
    .attr('class', 'networkGraph')
    .attr('transform', transform);

  // if we don't have any nodes
  if (!nodes) {
    networkGroup.selectAll('circle.nodes').remove();

    networkGroup.selectAll('path.links').remove();
    return;
  }

  //filter nodes to be those viewable with the current transform
  const viewableNodes = nodes.filter((node) => {
    const currentTransform = d3.zoomTransform(
      selection.node(),
    );
    const nodePosition = currentTransform.apply([
      node.x,
      node.y,
    ]);

    // The node is visible if its position falls within the SVG view box
    return (
      0 <= nodePosition[0] &&
      nodePosition[0] <= selection.attr('width') &&
      0 <= nodePosition[1] &&
      nodePosition[1] <= selection.attr('height')
    );
  });

  const viewablePacks = {};
  viewableNodes.forEach((node) => {
    const thisPack = d3
      .pack()
      .size([node.r * 2, node.r * 2])(packData[node.id]);
    thisPack.x = node.x;
    thisPack.y = node.y;
    //thisPack.r = node.r;
    // thisPack.descendants().forEach((child) => {
    //   child.x = child.x + node.x;
    //   child.y = child.y + node.y;
    // });
    viewablePacks[node.id] = thisPack;
  });

  !transform || transform.k <= zoomNameLimit
    ? networkGroup
        .selectAll('g.nodeGroup')
        .data(viewableNodes, (d) => d.id)
        .join(
          (enter) => {
            const nodeGroup = enter
              .append('g')
              .attr('class', 'nodeGroup');

            nodeGroup.selectAll('circle').remove();
            nodeGroup
              .append('circle')
              .attr('class', 'nodes')
              .attr('cx', (d) => d.x)
              .attr('cy', (d) => d.y)
              .attr('r', (d) => d.r)
              .attr('fill', (d) => d.fill)
              .attr('stroke', (d) => d.stroke)
              .attr('opacity', (d) =>
                (hoveredCluster !== null &&
                  hoveredCluster !== undefined) ||
                (clickedCluster !== null &&
                  clickedCluster !== undefined)
                  ? d.cluster === hoveredCluster ||
                    d.cluster === clickedCluster
                    ? 1
                    : 0.2
                  : 1,
              )
              .on('click', (event, d) => {
                setClickedNode(d);
              });

            return nodeGroup;
          },
          (update) => {
            update.selectAll('circle.descendants').remove();

            return update
              .select('circle.nodes')
              .attr('cx', (d) => d.x)
              .attr('cy', (d) => d.y)
              .attr('r', (d) => d.r)
              .attr('fill', (d) => d.fill)
              .attr('opacity', (d) =>
                (hoveredCluster !== null &&
                  hoveredCluster !== undefined) ||
                (clickedCluster !== null &&
                  clickedCluster !== undefined)
                  ? d.cluster === hoveredCluster ||
                    d.cluster === clickedCluster
                    ? 1
                    : 0.2
                  : 1,
              );
          },
          (exit) => exit.remove(),
        )
        .attr('stroke-width', nodeStrokeWidth)
        //.call(dragBehavior)
        .on('mouseover', function (event, d) {
          const [pointerX, pointerY] = d3.pointer(event);
          hmToolTip
            .transition()
            .duration(400)
            .style('opacity', 0.9);
          hmToolTip
            .html(
              `ID: ${d.id} <br>
        Name: ${d.text}`,
            )
            .style('position', 'fixed')
            .style('left', event.pageX + 10 + 'px')
            .style('top', event.pageY + 10 + 'px')
            .style('background', '#DDDDDD');
        })
        .on('mouseout', function () {
          hmToolTip
            .transition()
            .duration(500)
            .style('opacity', 0);
          hmToolTip.html('');
        })
        .style('cursor', 'pointer')
    : networkGroup
        .selectAll('g.nodeGroup')
        .data(viewableNodes, (d) => d.id)
        .join(
          (enter) => {
            const nodeGroup = enter
              .append('g')
              .attr('class', 'nodeGroup');

            nodeGroup.selectAll('circle').remove();

            nodeGroup
              .append('circle')
              .attr('class', 'nodes')
              .attr('cx', (d) => viewablePacks[d.id].x)
              .attr('cy', (d) => viewablePacks[d.id].y)
              .attr('r', (d) => viewablePacks[d.id].r)
              .attr('fill', '#777777')
              .attr('stroke', (d) => d.stroke);

            //add the descendants of viewablePacks[d.id]
            nodeGroup
              .selectAll('circle.descendants')
              .data(
                (d) => viewablePacks[d.id].children,
                (d) => d.id,
              )
              .join('circle')
              .attr('class', 'descendants')
              .attr(
                'cx',
                (d) => d.parent.x + d.x - d.parent.r,
              )
              .attr(
                'cy',
                (d) => d.parent.y + d.y - d.parent.r,
              )
              .attr('r', (d) => d.r)
              .attr('fill', (d) => d.data.fill)
              .attr('stroke', (d) => d.stroke)
              .attr('opacity', (d) =>
                (hoveredCluster !== null &&
                  hoveredCluster !== undefined) ||
                (clickedCluster !== null &&
                  clickedCluster !== undefined)
                  ? d.parent.cluster === hoveredCluster ||
                    d.parent.cluster === clickedCluster
                    ? 1
                    : 0.2
                  : 1,
              );

            return nodeGroup;
          },
          (update) => {
            update
              .select('circle.nodes')
              .attr('cx', (d) => viewablePacks[d.id].x)
              .attr('cy', (d) => viewablePacks[d.id].y)
              .attr('r', (d) => viewablePacks[d.id].r)
              .attr('fill', '#777777')
              .attr('stroke', (d) => d.stroke)
              .attr('opacity', (d) =>
                (hoveredCluster !== null &&
                  hoveredCluster !== undefined) ||
                (clickedCluster !== null &&
                  clickedCluster !== undefined)
                  ? d.cluster === hoveredCluster ||
                    d.cluster === clickedCluster
                    ? 1
                    : 0.2
                  : 1,
              )
              .call(dragBehavior)
              .style('cursor', 'pointer');

            update
              .selectAll('circle.descendants')
              .data(
                (d) => viewablePacks[d.id].children,
                (d) => d.id,
              )
              .join('circle')
              .attr('class', 'descendants')
              .attr(
                'cx',
                (d) => d.parent.x + d.x - d.parent.r,
              )
              .attr(
                'cy',
                (d) => d.parent.y + d.y - d.parent.r,
              )
              .attr('r', (d) => d.r)
              .attr('fill', (d) => d.data.fill)
              .attr('stroke', (d) => d.stroke)
              .attr('opacity', (d) =>
                (hoveredCluster !== null &&
                  hoveredCluster !== undefined) ||
                (clickedCluster !== null &&
                  clickedCluster !== undefined)
                  ? d.cluster === hoveredCluster ||
                    d.cluster === clickedCluster
                    ? 1
                    : 0.2
                  : 1,
              );

            return update;
          },
          (exit) => exit.remove(),
        )
        .attr('stroke-width', nodeStrokeWidth)
        .style('cursor', 'pointer')
        .call(dragBehavior)
        .on('click', (event, d) => {
          //event.originalEvent.stopPropogation();
          setClickedNode(d);
        })
        .on('mouseover', function (event, d) {
          const [pointerX, pointerY] = d3.pointer(event);
          hmToolTip
            .transition()
            .duration(200)
            .style('opacity', 0.9);
          hmToolTip
            .html(
              `ID: ${d.id} <br>
        Name: ${d.text}`,
            )
            .style('position', 'fixed')
            .style('left', event.pageX + 10 + 'px')
            .style('top', event.pageY + 10 + 'px')
            .style('background', '#DDDDDD');
        })
        .on('mouseout', function () {
          hmToolTip
            .transition()
            .duration(500)
            .style('opacity', 0);
          hmToolTip.html('');
        });

  // Define the arrowhead marker
  networkGroup
    .selectAll('marker.arrowhead')
    .data([null])
    .join('marker')
    .attr('id', 'arrowhead')
    .attr('class', 'arrowhead')
    .attr('markerWidth', 10)
    .attr('markerHeight', 10)
    .attr('refX', 9)
    .attr('refY', 5)
    .attr('orient', 'auto-start-reverse')
    .append('path')
    .attr('d', 'M0,1 L0,9 L10,5 z')
    .attr('fill', 'context-stroke'); // Path for arrowhead shape

  //thank ChatGPT for the line ends code
  networkGroup
    .selectAll('path.links')
    .data(links)
    .join('path')
    .attr('class', 'links')
    .attr('d', (d) => {
      if (d.source.id !== d.target.id) {
        const source = {
          x:
            d.source.x +
            (d.source.r + 1) *
              Math.cos(
                Math.atan2(
                  d.target.y - d.source.y,
                  d.target.x - d.source.x,
                ),
              ),
          y:
            d.source.y +
            (d.source.r + 1) *
              Math.sin(
                Math.atan2(
                  d.target.y - d.source.y,
                  d.target.x - d.source.x,
                ),
              ),
        };
        const target = {
          x:
            d.target.x -
            (d.target.r + 1) *
              Math.cos(
                Math.atan2(
                  d.target.y - d.source.y,
                  d.target.x - d.source.x,
                ),
              ),
          y:
            d.target.y -
            (d.target.r + 1) *
              Math.sin(
                Math.atan2(
                  d.target.y - d.source.y,
                  d.target.x - d.source.x,
                ),
              ),
        };
        const dx = target.x - source.x;
        const dy = target.y - source.y;
        const angle = Math.atan2(dy, dx); // Angle of the line segment
        const controlPoint = {
          x:
            source.x +
            dx / 2 +
            5 * Math.cos(angle + Math.PI / 2),
          y:
            source.y +
            dy / 2 +
            5 * Math.sin(angle + Math.PI / 2),
        }; // Control point 10 units above the line
        //return `M${source.x},${source.y} Q${controlPoint.x},${controlPoint.y} ${target.x},${target.y}`; // Quadratic Bézier curve path
        return `M${source.x},${source.y} ${target.x},${target.y}`; // straight path
      } else {
        return `M${d.source.x + d.source.r},${d.source.y} 
                Q${d.source.x + d.source.r + 15},${d.source.y + 11} , 
                ${d.source.x + d.source.r + 15},${d.source.y - 3}T 
                ${d.source.x + d.source.r},${d.source.y}`;
      }
    })
    .attr('fill', 'none')
    .attr('stroke', (d) => d.color)
    .attr('stroke-width', linkStrokeWidth)
    //.attr('marker-end', 'url(#arrowhead)')
    .attr('opacity', (d) =>
      (hoveredCluster !== null &&
        hoveredCluster !== undefined) ||
      (clickedCluster !== null &&
        clickedCluster !== undefined)
        ? d.source.cluster === hoveredCluster ||
          d.source.cluster === clickedCluster
          ? 1
          : 0.2
        : 1,
    );

  //transition code is based on https://vizhub.com/curran/parallel-coordinates-with-brushing?edit=files&file=parallelCoordinates.js
  const t = d3
    .transition()
    .duration(750)
    .ease(d3.easeLinear);

  !transform || transform.k <= zoomNameLimit
    ? networkGroup
        .selectAll('text.nodes')
        .data([])
        .join(
          (enter) => enter,
          (update) => update,
          (exit) =>
            exit.call((exit) =>
              exit
                .transition(t)
                .attr('opacity', 0)
                .remove(),
            ),
        )
    : networkGroup
        .selectAll('text.nodes')
        .data(viewableNodes, (d) => d.id)
        .join(
          (enter) =>
            enter
              .append('text')
              .attr('class', 'nodes')
              .attr('opacity', 0)
              .call((enter) =>
                enter.transition(t).attr('opacity', 1),
              ),
          (update) =>
            update.call((update) =>
              update.transition(t).attr('opacity', 1),
            ),
          (exit) => exit.remove(),
        )
        .attr('class', 'nodes')
        .attr('x', (d) => d.x - d.r / Math.sqrt(2) - 1)
        .attr('y', (d) => d.y + d.r / Math.sqrt(2) + 1)
        .attr('text-anchor', 'end')
        .attr('dominant-baseline', 'hanging')
        .attr('font-size', textFontSize / transform.k)
        .text((d) => d.text);
};
