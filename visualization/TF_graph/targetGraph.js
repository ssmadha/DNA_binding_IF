//import * as d3 from 'd3';

export const targetGraph = (
  selection,
  { nodes, clickedNode, links },
) => {
  if (!nodes || nodes.length === 0 || !links) {
    return;
  }

  const width = 250;
  const height = 250;
  const padding = 40;

  // Filter links to show only connections involving clicked nodes
  const clickedNodeId = clickedNode.map((d) => d.id);
  const miniGraphLinks = links.filter(
    (link) =>
      clickedNodeId.includes(link.source.id || link.source),
  );

  // Get all nodes involved in the mini graph (clicked nodes + connected nodes)
  const connectedNodeIds = new Set(clickedNode);
  miniGraphLinks.forEach((link) => {
    connectedNodeIds.add({
      id: link.target,
      text: link.text,
    });
  });

  // Find the actual node objects for all involved nodes
  // Create independent copies to avoid modifying the main graph nodes
  // Separate clicked nodes from connected nodes
  const clickedNodesMap = {};
  const connectedNodes = [];

  Array.from(connectedNodeIds).forEach((d) => {
    const existingNode = nodes.find((n) => n.id === d.id);
    const node = existingNode
      ? {
          ...existingNode,
          x: 0,
          y: 0,
        }
      : {
          id: d.id,
          text: d.text,
          x: 0,
          y: 0,
          r: 5,
          fill: '#999999',
          stroke: '#333333',
        };

    if (clickedNodeId.includes(d.id)) {
      clickedNodesMap[d.id] = node;
    } else {
      connectedNodes.push(node);
    }
  });

  const allNodes = [...Object.values(clickedNodesMap), ...connectedNodes];

  // console.log('Mini graph nodes:', allNodes);
  // console.log('Mini graph links:', miniGraphLinks);

  // Normalize link objects to have source/target as objects with id property
  const normalizedLinks = miniGraphLinks.map((link) => ({
    ...link,
    source: typeof link.source === 'string' ? { id: link.source } : link.source,
    target: typeof link.target === 'string' ? { id: link.target } : link.target,
  }));

  // Create or update SVG
  const svg = selection
    .selectAll('svg.targetGraphSvg')
    .data([null])
    .join('svg')
    .attr('class', 'targetGraphSvg')
    .attr('width', width)
    .attr('height', height)
    .style('border', '1px solid #ccc')
    .style('background', '#ffffff');

  // Create group for graph elements
  const graphGroup = svg
    .selectAll('g.graphGroup')
    .data([null])
    .join('g')
    .attr('class', 'graphGroup')
    .attr('transform', `translate(${padding},${padding})`);

  const graphWidth = width - 2 * padding;
  const graphHeight = height - 2 * padding;
  const centerX = graphWidth / 2;
  const centerY = graphHeight / 2;
  const radius = 80;

  // Position clicked node(s) at center
  Object.values(clickedNodesMap).forEach((node) => {
    node.x = centerX;
    node.y = centerY;
  });

  // Position connected nodes evenly around the center
  const numConnectedNodes = connectedNodes.length;
  if (numConnectedNodes > 0) {
    connectedNodes.forEach((node, i) => {
      const angle = (2 * Math.PI * i) / numConnectedNodes;
      node.x = centerX + radius * Math.cos(angle);
      node.y = centerY + radius * Math.sin(angle);
    });
  }

  // Constrain nodes to stay within bounds
  allNodes.forEach((node) => {
    const nodeRadius = node.r || 5;
    node.x = Math.max(nodeRadius, Math.min(graphWidth - nodeRadius, node.x));
    node.y = Math.max(nodeRadius, Math.min(graphHeight - nodeRadius, node.y));
  });

  // Draw links
  const linkSelection = graphGroup
    .selectAll('line.miniLink')
    .data(normalizedLinks, (d) => `${d.source.id}-${d.target.id}`)
    .join('line')
    .attr('class', 'miniLink')
    .attr('x1', (d) => allNodes.find((n) => n.id === d.source.id).x)
    .attr('y1', (d) => allNodes.find((n) => n.id === d.source.id).y)
    .attr('x2', (d) => allNodes.find((n) => n.id === d.target.id).x)
    .attr('y2', (d) => allNodes.find((n) => n.id === d.target.id).y)
    .attr('stroke', (d) => d.color || '#7e7e7e')
    .attr('stroke-width', 1.4)
    .attr('opacity', 0.6);

  // Draw nodes
  const nodeSelection = graphGroup
    .selectAll('circle.miniNode')
    .data(allNodes, (d) => d.id)
    .join('circle')
    .attr('class', 'miniNode')
    .attr('cx', (d) => d.x)
    .attr('cy', (d) => d.y)
    .attr('r', (d) => d.r)
    .attr('fill', (d) => d.fill)
    .attr('stroke', (d) => d.stroke)
    .attr('stroke-width', (d) =>
      clickedNodeId.includes(d.id) ? 2 : 1,
    );

  // Draw labels
  const labelSelection = graphGroup
    .selectAll('text.miniNodeLabel')
    .data(allNodes, (d) => d.id)
    .join('text')
    .attr('class', 'miniNodeLabel')
    .attr('x', (d) => d.x + d.r/Math.sqrt(2) + 4)
    .attr('y', (d) => d.y + d.r/Math.sqrt(2) + 8)
    .attr('font-size', 10)
    .attr('fill', '#333333')
    .text((d) => d.text);

  // Add tooltip
  let tooltip = d3.select('body .tooltip-mini');
  if (tooltip.empty()) {
    tooltip = d3
      .select('body')
      .append('div')
      .attr('class', 'tooltip-mini')
      .style('opacity', 0);
  }

  nodeSelection
    .on('mouseover', function (event, d) {
      tooltip
        .transition()
        .duration(200)
        .style('opacity', 0.9);
      tooltip
        .html(`ID: ${d.id}<br/>Name: ${d.text}`)
        .style('position', 'fixed')
        .style('left', event.pageX + 10 + 'px')
        .style('top', event.pageY + 10 + 'px')
        .style('background', '#DDDDDD')
        .style('padding', '8px')
        .style('border-radius', '4px')
        .style('font-size', '12px');
    })
    .on('mouseout', function () {
      tooltip
        .transition()
        .duration(500)
        .style('opacity', 0);
    });
};
