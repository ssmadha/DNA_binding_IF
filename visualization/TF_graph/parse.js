//import * as d3 from 'd3';

export const parseData = ({
  data,
  id,
  text,
  simMinX,
  simMaxX,
  simMinY,
  simMaxY,
  radiusValue,
  radiusExtent,
  nodeFillValue,
  nodeStrokeValue,
}) => {
  const radiusScale = d3.scaleLinear(
    d3.extent(data, radiusValue),
    radiusExtent,
  );

  const colorFillScale = d3.scaleSequential(
    [0, d3.max(data, nodeFillValue)],
    d3.interpolateInferno,
  );

  const colorStrokeScale = d3.scaleOrdinal(
    Object.values(nodeStrokeValue),
    d3.schemePaired,
  );

  const xRanks = d3.rank(
    data,
    (a, b) =>
      d3.ascending(nodeFillValue(a), nodeFillValue(b)) || 1,
  );
  const yRanks = d3.rank(
    data,
    (a, b) =>
      d3.ascending(radiusValue(a), radiusValue(b)) || 1,
  );

  console.log(xRanks);
  return {
    nodes: data.map((d, i) => {
      return {
        id: id !== undefined ? id(d) : i,
        x:
          (xRanks[i] / data.length) * (simMaxX - simMinX) +
          simMinX,
        y:
          (yRanks[i] / data.length) * (simMaxY - simMinY) +
          simMinY,
        text: text(d),
        r: radiusScale(radiusValue(d)),
        fill: colorFillScale(nodeFillValue(d)),
        stroke: colorStrokeScale(nodeStrokeValue[id(d)]),
        cluster: nodeStrokeValue[id(d)],
      };
    }),
    radiusScale,
    colorFillScale,
    colorStrokeScale,
  };
};

export const parsePackData = ({
  data,
  dataID,
  subData,
  parentID,
  childID,
  radiusValue,
  colorValue,
}) => {
  const colorScale = d3.scaleSequential(
    [0, d3.max(data, colorValue)],
    d3.interpolateInferno,
  );
  const radiusScale = d3.scaleLinear();
  let result = {};
  //iterate over the values in the "Gene" column of data
  data.map((d) => {
    const subDataFiltered = subData.filter(
      (d1) => parentID(d1) === dataID(d),
    );

    subDataFiltered.forEach((d1) => {
      d1.fill = colorScale(colorValue(d1));
    });
    //append a value to subDataFiltered where parentID = "" and childID = root where parentID and childID are accessor functions
    subDataFiltered.unshift({
      parent: '',
      id: dataID(d),
    });

    const root = d3
      .stratify()
      .id((d1, i) => (i === 0 ? d1.id : childID(d1)))
      .parentId((d1, i) =>
        i === 0 ? d1.parent : parentID(d1),
      )(subDataFiltered);

    root.sum((d1) =>
      radiusValue(d1) ? radiusValue(d1) : 0,
    );
    result[dataID(d)] = root;
  });
  return result;
};

export const parseEdges = ({
  data,
  dataID,
  edges,
  sourceValue,
  targetValue,
  edgeColorValue,
}) => {
  const colorScale = d3.scaleLinear(
    d3.extent(edges, edgeColorValue),
    ['grey', 'red'],
  );
  return edges
    .map((d) => {
      // check if the source and target of d are in data's id column
      const sourceExists = data.some(
        (item) => dataID(item) === sourceValue(d),
      );
      const targetExists = data.some(
        (item) => dataID(item) === targetValue(d),
      );
      if (sourceExists && targetExists) {
        return {
          source: sourceValue(d),
          target: targetValue(d),
          color: colorScale(edgeColorValue(d)),
        };
      } else {
        return null;
      }
    })
    .filter((d) => d !== null);
};
