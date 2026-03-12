//import * as d3 from '../d3';

import { graph } from './graph.js';
import { body } from './body.js';
 import { bodyData, organDict } from '../body/index.js';
import {
  asif_melted,
  gene_cluster,
  hi_union,
  transcript_data,
} from '../asif/index.js';
import {
  parseData,
  parsePackData,
  parseEdges,
} from './parse.js';
import { simulate } from './simulate.js';
import { clusterLegend } from './clusterLegend.js';
import { nodeLegend } from './nodeLegend.js';
import { targetGraph } from './targetGraph.js';
import { infoView } from './infoView.js';


export const main = (container, { state, setState }) => {
  // compute width/height only once and store in state to avoid the growth
  // feedback loop described above. if the container does not have a fixed
  // height, `container.clientHeight` will grow whenever we draw an svg with
  // that same height, so subsequent renders keep inflating the page.
  let width = state.dimensions?.width;
  let height = state.dimensions?.height;

  if (width === undefined || height === undefined) {
    width = container.clientWidth;
    // fall back to viewport size if the container is still empty
    height = container.clientHeight || window.innerHeight * 0.8;

    setState((s) => ({
      ...s,
      dimensions: { width, height },
    }));
  }

  const {
    transform,
    nodes,
    links,
    radiusScale,
    colorFillScale,
    colorStrokeScale,
    packData,
    keepData,
    freezeSim,
    clickedCluster,
    hoveredCluster,
    clickedOrgan,
    hoveredOrgan,
    clickedNodes: clickedNode,
    filterLimit,
    geneSearch,
  } = state;

  if (filterLimit === undefined) {
    setState((state) => ({
      ...state,
      filterLimit: 20,
    }));
    return;
  }

  const nodeLegendWidth = 100;
  const nodeLegendHeight = 500;

  const bodyWidth = 360;
  const bodyHeight = 550;

  const netGraphX = bodyWidth;
  const netGraphY = 0;
  const netGraphWidth = width - netGraphX - nodeLegendWidth;
  const netGraphHeight = bodyHeight;

  const simMinX = 10;
  const simMaxX = netGraphWidth - 10;
  const simMinY = 10;
  const simMaxY = netGraphHeight - 10;

  //zooming in and out once gets the transform scale to 0.999...
  const zoomNameLimit = 1.0;

  //based on https://vizhub.com/curran/multiple-menus-with-d3 and https://vizhub.com/curran/menu-with-d3
  const containerSelection = d3
    .select(container)
    .attr('class', 'ui-elements-container');

  const menusContainer = containerSelection
    .selectAll('div.menus-container')
    .data([null])
    .join('div')
    .attr('class', 'menus-container row g-6');

  const filterRangeDiv = menusContainer
    .selectAll('div.filter-range-container')
    .data([null])
    .join('div')
    .attr('class', 'filter-range-container col');

  filterRangeDiv
    .selectAll('label')
    .data([null])
    .join('label')
    .attr('class', 'form-label')
    .attr('for', 'filterRange')
    .text(`Expression Limit: ${filterLimit}`);

  filterRangeDiv
    .selectAll('input')
    .data([null])
    .join('input')
    .attr('type', 'range')
    .attr('class', 'form-range')
    .attr('id', 'filterRange')
    .attr('min', 0)
    .attr(
      'max',
      nodes
        ? d3.max(nodes, (d) => d.Expression)
        : d3.max(asif_melted, (d) => d.Expression),
    )
    .attr('value', filterLimit)
    .on('input', (event) =>
      setState((state) => ({
        ...state,
        filterLimit: event.target.value,
        freezeSim: true,
      })),
    )
    .on('change', (event) =>
      setState((state) => ({
        ...state,
        filterLimit: event.target.value,
        keepData: false,
        freezeSim: false,
      })),
    );

  const geneSearchDiv = menusContainer
    .selectAll('div.gene-search-container')
    .data([null])
    .join('div')
    .attr('class', 'gene-search-container col');

  geneSearchDiv
    .selectAll('label')
    .data([null])
    .join('label')
    .attr('class', 'form-label')
    .attr('for', 'geneSearch')
    .text(`Gene Search:`);

  geneSearchDiv
    .selectAll('input')
    .data([null])
    .join('input')
    .attr('type', 'text')
    .attr('class', 'form-text')
    .attr('id', 'geneSearch')
    .on('input', (event) => {
      const searchText = event.target.value;
      console.log(searchText);
      if (searchText || searchText != '') {
        setState((state) => ({
          ...state,
          geneSearch: searchText,
        }));
      } else {
        setState((state) => ({
          ...state,
          geneSearch: undefined,
        }));
      }
    });

  const setClickedCluster = (clickedCluster) => {
    setState((state) => ({
      ...state,
      clickedCluster,
      keepData: true,
      freezeSim: true,
    }));
  };

  const setHoveredCluster = (hoveredCluster) => {
    setState((state) => ({
      ...state,
      hoveredCluster,
      keepData: true,
      freezeSim: true,
    }));
  };

  const clusterLegendDiv = menusContainer
    .selectAll('div.cluster-legend')
    .data([null])
    .join('div')
    .attr('class', 'cluster-legend col');

  if (nodes) {
    const clusterLegendSvg = clusterLegendDiv
      .selectAll('svg.cluster-legend')
      .data([null])
      .join('svg')
      .attr('class', 'cluster-legend')
      .attr('width', 360)
      .attr('height', 50);

    clusterLegend(clusterLegendSvg, {
      colorScale: colorStrokeScale,
      legendLabel: 'Cluster',
      legendLabelY: 20,
      legendPadding: 70,
      tickSpacing: 25,
      legendCircleRadius: 5,
      tickPadding: 10,
      hoveredCluster,
      setHoveredCluster,
      clickedCluster,
      setClickedCluster,
    });
  }

  const setClickedOrgan = (clickedOrgan) => {
    setState((state) => ({
      ...state,
      transform: undefined,
      clickedOrgan,
      keepData: false,
      freezeSim: false,
    }));
  };

  const setHoveredOrgan = (hoveredOrgan) => {
    setState((state) => ({
      ...state,
      hoveredOrgan,
      keepData: true,
      freezeSim: true,
    }));
  };

  const viewDiv = containerSelection
    .selectAll('div.view-container')
    .data([null])
    .join('div')
    .attr('class', 'flex-row view-container');

  const bodyGroup = viewDiv
    .selectAll('div.body')
    .data([null])
    .join('div')
    .attr('class', 'body');

  const bodySVG = bodyGroup
    .selectAll('svg.body')
    .data([null])
    .join('svg')
    .attr('class', 'body')
    .attr('width', bodyWidth)
    .attr('height', bodyHeight)
    //.attr('preserveAspectRatio', 'none')
    .style('background', null);

  bodySVG.call(body, {
    bodyData,
    organDict,
    legendLabel: 'Organ',
    legendLabelFontSize: 22,
    legendX: 20,
    legendY: 25,
    legendPadding: 13,
    legendFontSize: 12,
    tickSpacing: 13,
    tickPadding: 0,
    clickedOrgan,
    setClickedOrgan,
    hoveredOrgan,
    setHoveredOrgan,
  });

  const tissue = clickedOrgan;

  if (!keepData) {
    if (clickedOrgan) {
      const filteredData = asif_melted.filter(
        (d) =>
          d.Tissue === tissue &&
          d.Expression >= filterLimit,
      );
      const {
        nodes,
        radiusScale,
        colorFillScale,
        colorStrokeScale,
      } = parseData({
        data: filteredData,
        id: (d) => d['Gene ID'],
        text: (d) => d['Gene Name'],
        simMinX,
        simMaxX,
        simMinY,
        simMaxY,
        radiusValue: (d) => d.Expression,
        radiusExtent: [5, 15],
        nodeFillValue: (d) => d.ASIF,
        nodeStrokeValue: gene_cluster[0],
      });
      setState((state) => ({
        ...state,
        nodes,
        radiusScale,
        colorFillScale,
        colorStrokeScale,
        links: parseEdges({
          data: filteredData,
          dataID: (d) => d['Gene ID'],
          edges: hi_union,
          sourceValue: (d) => d.source,
          targetValue: (d) => d.target,
          edgeColorValue: (d) => d.coverage,
        }),
        packData: parsePackData({
          data: filteredData,
          dataID: (d) => d['Gene ID'],
          subData: transcript_data.filter(
            (d) => d.Tissue == tissue,
          ),
          parentID: (d) => d.Gene,
          childID: (d) => d.Transcript,
          radiusValue: (d) => d.Expression,
          colorValue: (d) => d.ASIF,
        }),
        keepData: true,
      }));
      return;
    } else {
      setState((state) => ({
        ...state,
        nodes: null,
        links: null,
        keepData: true,
      }));
    }
  }

  const graphDiv = viewDiv
    .selectAll('div.graph-container')
    .data([null])
    .join('div')
    .attr('class', 'graph-container');

  const graphSvg = graphDiv
    .selectAll('svg.graph')
    .data([null])
    .join('svg')
    .attr('class', 'graph')
    .attr('x', netGraphX)
    .attr('y', netGraphY)
    .attr('width', netGraphWidth)
    .attr('height', netGraphHeight);

  // filter nodes based on those that exist in links
  const filteredNodes = nodes
    ? nodes.filter((node) =>
        links.some(
          (link) =>
            link.source === node.id ||
            link.target === node.id ||
            link.source.id === node.id ||
            link.target.id === node.id,
        ),
      )
    : [];

  const filteredLinks = true ? links : [];

  if (!freezeSim && clickedOrgan) {
    simulate({
      nodes: filteredNodes,
      links: filteredLinks,
      netGraphWidth,
      netGraphHeight,
      simMinX,
      simMaxX,
      simMinY,
      simMaxY,
      linkDistance: 30,
      linkStrength: 0.9,
      linkIterations: 10,
      manyBodyStrength: -150,
      manyBodyMinDist: 1,
      manyBodyMaxDist: (d) => d.r + 60,
      collideRadius: (d) => d.r + 15,
      centerX: netGraphWidth / 2,
      centerY: netGraphHeight / 2,
      centerStrength: 0.6,
      strengthXY: 0.5,
      velocityDecay: 0.8,
    });
  }
  console.log(links);
  //graphSvg.attr('transform', transform);

  //zoom code from https://vizhub.com/curran/d6f1170765c84a498caa6ea11403e3be
  const zoomBehavior = d3
    .zoom()
    .scaleExtent([0.75, 15])
    .on('zoom', (event) => {
      setState((state) => ({
        ...state,
        transform: event.transform,
        freezeSim: true,
      }));
    });

  const setClickedNode = (clickedNode) => {
    setState((state) => ({
      ...state,
      clickedNodes: [clickedNode],
      keepData: true,
      freezeSim: true,
    }));
    return;
    if (!clickedNode || clickedNode.length == 0) {
      // console.log(clickedNodes);
      // console.log('empty clicked nodes');
      setState((state) => ({
        ...state,
        clickedNodes: [clickedNode],
        keepData: true,
        freezeSim: true,
      }));
    } else if (clickedNode.includes(clickedNode)) {
      // console.log('removing');
      setState((state) => ({
        ...state,
        clickedNodes: clickedNode.filter(
          (d) => d.id !== clickedNode.id,
        ),
        keepData: true,
        freezeSim: true,
      }));
    } else {
      // console.log('pushing');
      clickedNode.push(clickedNode);
      setState((state) => ({
        ...state,
        clickedNodes: [clickedNode],
        keepData: true,
        freezeSim: true,
      }));
    }
  };

  graphSvg.call(zoomBehavior);

  //console.log(nodes);
  graph(graphSvg, {
    nodes: geneSearch
      ? nodes.filter((d) =>
          d.text
            .toLowerCase()
            .includes(geneSearch.toLowerCase()),
        )
      : filteredNodes,
    links: filteredLinks,
    packData,
    nodeStrokeWidth: 3,
    linkStrokeWidth: 1,
    textFontSize: 18,
    transform,
    zoomNameLimit,
    hoveredCluster,
    clickedCluster,
    setClickedNode,
    setState,
  });

  const nodeLegendDiv = viewDiv
    .selectAll('div.nodeLegend')
    .data([null])
    .join('div')
    .attr('class', 'nodeLegend');

  const nodeLegendSvg = nodeLegendDiv
    .selectAll('svg.nodeLegend')
    .data([null])
    .join('svg')
    .attr('class', 'nodeLegend')
    .attr('width', nodeLegendWidth)
    .attr('height', nodeLegendHeight);

  if (nodes) {
    nodeLegend(nodeLegendSvg, {
      colorScale: colorFillScale,
      colorLabel: 'ASIF score',
      sizeScale: radiusScale,
      sizeLabel: 'Expression',
      legendLabelFontSize: 16,
      legendPadding: 20,
      legendX: 10,
      legendY: 15,
      legendFontSize: 12,
      colorNumber: 100,
      colorPadding: 15,
      colorHeight: 250,
      colorWidth: 20,
      colorTickNumber: 8,
      legendSpacing: 35,
      sizeNumber: 4,
      sizeSpacing: 30,
      sizePadding: 17,
    });
  }


  const clickedInfoDiv = containerSelection
    .selectAll('div.clicked-info-container')
    .data([null])
    .join('div')
    .attr('class', 'flex-row clicked-info-container');

  const targetDiv = clickedInfoDiv
    .selectAll('div.targetGraph')
    .data([null])
    .join('div')
    .attr('class', 'targetGraph');

  const infoDiv = clickedInfoDiv
    .selectAll('div.info-container')
    .data([null])
    .join('div')
    .attr('class', 'info-container');

  if (clickedNode) {
    targetGraph(targetDiv, {
      nodes,
      clickedNode,
      links,
      data: asif_melted,
      subData: transcript_data,
      tissue,
    });

    infoView(infoDiv, {
      nodes: clickedNode,
      data: asif_melted,
      subData: transcript_data,
      tissue,
    });
  }
};

const container = document.getElementById("viz");

let state = {};

const setState = (updater) => {
  state = updater(state);
  main(container, { state, setState });
};

main(container, {
  state,
  setState
});