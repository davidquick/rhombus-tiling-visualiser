class Point {
    constructor (x, y) {
        this.x = x;
        this.y = y;
    }

    toString() {
        return `Point(${this.x}, ${this.y})`
    }
}

class Vertex {
    constructor(n, previous_bearing) {
        this.n = n;
        this.previous_bearing = previous_bearing;
    }

    get next_bearing() {
        return this.previous_bearing + 1
    }

    get neighbouring_bearings() {
        return [this.previous_bearing, this.next_bearing]
    }

    get x() {
        return Math.cos(Math.PI * (this.previous_bearing+0.5)/this.n)
    }

    get y() {
        return Math.sin(Math.PI * (this.previous_bearing+0.5)/this.n)
    }

    get point() {
        return new Point(this.x, this.y)
    }

    get key() {
        return [this.previous_bearing % this.n, this.previous_bearing < this.n]
    }

    toString() {
        return `Vertex(${this.previous_bearing}.5/${this.n}, ${this.x}, ${this.y})`
    }
}

function orientation_colour(orientation, n, colouring_type, colour_base) {
    var base_hsl = d3.hsl(colour_base);
    var value = orientation/n
    switch (colouring_type) {
        case "hue 360":
            return d3.rgb(d3.hsl(base_hsl.h + (360*value), base_hsl.s, base_hsl.l));
        case "hue 360 symmetric":
            return d3.rgb(d3.hsl(base_hsl.h + (720*(value < 0.5 ? value : 1 - value)), base_hsl.s, base_hsl.l));
        case "hue 90":
            return d3.rgb(d3.hsl(base_hsl.h + (90*value), base_hsl.s, base_hsl.l));
        case "hue 90 symmetric":
            return d3.rgb(d3.hsl(base_hsl.h + (180*(value < 0.5 ? value : 1 - value)), base_hsl.s, base_hsl.l));
        case "saturation":
            return d3.rgb(d3.hsl(base_hsl.h, value, base_hsl.l));
        case "luminance":
            return d3.rgb(d3.hsl(base_hsl.h, base_hsl.s, value));
        case "constant":
        default:
            return base_hsl;
    }
}

function orientations_colour(orientation_1, orientation_2, n, edge_colour_type, edge_colour_base, face_colour_type, face_colour_base) {
    var edge_base_hsl = d3.hsl(edge_colour_base);
    var face_base_hsl = d3.hsl(face_colour_base);

    var average_orientation = (orientation_1 + orientation_2)/2;

    switch (face_colour_type) {
        case "colour blend":
            var colour_1 = d3.rgb(orientation_colour(orientation_1, n, edge_colour_type, edge_colour_base));
            var colour_2 = d3.rgb(orientation_colour(orientation_2, n, edge_colour_type, edge_colour_base));

            return d3.rgb(Math.sqrt(colour_1.r**2 + colour_2.r**2), Math.sqrt(colour_1.g**2 + colour_2.g**2), Math.sqrt(colour_1.b**2 + colour_2.b**2));
        case "average cyclic":
            if (Math.abs(average_orientation + (this.n/2) - orientation_2) % this.n < Math.abs(average_orientation - orientation_2) % this.n) {
                average_orientation = average_orientation + (this.n/2);
            }
        case "average edge":
            return orientation_colour(average_orientation, n, edge_colour_type, edge_colour_base);
        case "constant":
        default:
            return face_base_hsl;
    }
}

class Line {
    constructor (n, start, end, orientation, step, placed) {
        this.n = n;
        this.start = start;
        this.end = end;
        this.orientation = orientation;
        this.step = step;
        this.placed = placed;
    }

    get key() {
        return [this.orientation, this.step, this.placed];
    }

    colour(edge_colouring, edge_colour_base) {
        return orientation_colour(this.orientation, this.n, edge_colouring, edge_colour_base)
    }

    get line() {
        return this;
    }

    toString() {
        return `Line(${this.orientation}/${this.n}, ${this.start.toString()}, ${this.end.toString()})`
    }

    x1(x_scale) {
        var x1 = this.start.x;
        if (this.placed) {
            return x_scale(x1);
        } else {
            var scale_factor = 4
            return x_scale((x1/scale_factor - Math.min(this.step, this.orientation)/this.n));
        }
    }

    x1_collapsed(x_scale) {
        var x1 = (this.start.x + this.end.x)/2;
        if (this.placed) {
            return x_scale(x1);
        } else {
            var scale_factor = 4
            return x_scale((x1/scale_factor - Math.min(this.step, this.orientation)/this.n));
        }
    }

    x2(x_scale) {
        var x2 = this.end.x;
        if (this.placed) {
            return x_scale(x2);
        } else {
            var scale_factor = 4
            return x_scale((x2/scale_factor - Math.min(this.step, this.orientation)/this.n));
        }
    }

    x2_collapsed(x_scale) {
        var x2 = (this.start.x + this.end.x)/2;
        if (this.placed) {
            return x_scale(x2);
        } else {
            var scale_factor = 4
            return x_scale((x2/scale_factor - Math.min(this.step, this.orientation)/this.n));
        }
    }

    y1(y_scale) {
        var y1 = this.start.y;
        if (this.placed) {
            return y_scale(y1);
        } else {
            var scale_factor = 4
            return y_scale(1 + (y1/scale_factor - Math.max(this.step, this.orientation)/this.n));
        }
    }

    y1_collapsed(y_scale) {
        var y1 = (this.start.y + this.end.y)/2;
        if (this.placed) {
            return y_scale(y1);
        } else {
            var scale_factor = 4
            return y_scale(1 + (y1/scale_factor - Math.max(this.step, this.orientation)/this.n));
        }
    }

    y2(y_scale) {
        var y2 = this.end.y;
        if (this.placed) {
            return y_scale(y2);
        } else {
            var scale_factor = 4
            return y_scale(1 + (y2/scale_factor - Math.max(this.step, this.orientation)/this.n));
        }
    }

    y2_collapsed(y_scale) {
        var y2 = (this.start.y + this.end.y)/2;
        if (this.placed) {
            return y_scale(y2);
        } else {
            var scale_factor = 4
            return y_scale(1 + (y2/scale_factor - Math.max(this.step, this.orientation)/this.n));
        }
    }
}

class Edge {
    constructor (n, bearing) {
        this.n = n;
        this.bearing = bearing;
    }

    get orientation() {
        return this.bearing % n;
    }

    get start() {
        return new Vertex(this.n, this.bearing - 1);
    }

    get end() {
        return new Vertex(this.n, this.bearing);
    }

    get line() {
        return new Line(this.n, this.start, this.end, this.bearing % this.n, 0, true)
    }

    colour(edge_colouring, edge_colour_base) {
        return orientation_colour(this.orientation, this.n, edge_colouring, edge_colour_base)
    }

    get key() {
        return [this.bearing % this.n, this.bearing < this.n]
    }

    get placed() {
        return true;
    }

    toString() {
        return `Edge(${this.bearing}/${this.n}, ${this.start.toString()}, ${this.end.toString()})`
    }
}

class Tile {
constructor(n, orientation_1, orientation_2) {
        this.n = n;
        this.orientation_1 = Math.min(orientation_1, orientation_2);
        this.orientation_2 = Math.max(orientation_1, orientation_2);
    }

    toString() {
        return `Tile(${this.n}, ${this.orientation_1}x${this.orientation_2})`
    }
}

class Placement {
    constructor(n, bearing_1, bearing_2) {
        this.n = n;
        this.bearing_1 = bearing_1;
        this.bearing_2 = bearing_2;
    }

    get tile() {
        return new Tile(this.n, this.bearing_1 % this.n, this.bearing_2 % this.n)
    }

    toString() {
        return `Placement(${this.n}, ${this.bearing_1}x${this.bearing_2})`
    }
}

class Face {
    constructor (n, placement, points, placed, correct_order=true) {
        this.n = n;
        this.placement = placement;
        this.points = points;
        this.placed = placed == null ? true : false;

        if (correct_order) {
            var ordered_bearings = [placement.bearing_1, placement.bearing_2, (placement.bearing_1 + this.n) % (2*this.n), (placement.bearing_2 + this.n) % (2*this.n)];
            var min_bearing = ordered_bearings.reduce(function(lowest, next, index) { return next < ordered_bearings[lowest] ? index : lowest; }, 0);
            this.points = (this.points.concat(this.points)).splice(min_bearing, this.points.length);
        }
    }

    get key() {
        var tile = this.placement.tile;
        // console.log([tile.orientation_1, tile.orientation_2])
        return [tile.orientation_1, tile.orientation_2]
    }

    colour(edge_colour_type, edge_colour_base, face_colour_type, face_colour_base) {
        var tile = this.placement.tile;
        return orientations_colour(tile.orientation_1, tile.orientation_2, this.n, edge_colour_type, edge_colour_base, face_colour_type, face_colour_base)
    }

    get ordered_points() {
        return this.points;
    }

    shape(x_scale, y_scale) {
        var face = this;
        return "M " + this.ordered_points.map( function(p) {
            if (face.placed) {
                return [x_scale(p.x), y_scale(p.y)].join(",");
            } else {
                var scale_factor = 4
                return [x_scale((p.x/scale_factor - face.placement.bearing_1/face.n)), y_scale(1 + (p.y/scale_factor - face.placement.bearing_2/face.n))].join(",");
            }
        }).join(" L ") + " Z ";
    }

    tiny_shape(x_scale, y_scale) {
        var centre_x = this.ordered_points.map(p => p.x).reduce((a, b) => (a + b)) /4;
        var centre_y = this.ordered_points.map(p => p.y).reduce((a, b) => (a + b)) /4;
        var tiny_points = this.ordered_points.map(p => new Point(0.999*centre_x + 0.001*p.x, 0.999*centre_y + 0.001*p.y));
        return "M " + tiny_points.map( function(p) {
            return [x_scale(p.x), y_scale(p.y)].join(",");
        }).join(" L ") + " Z ";
    }

    toString() {
        return `Face(${this.n}, ${this.placement.bearing_1}x${this.placement.bearing_2})`
    }
}

class BearingCycle {
    constructor (n) {
        this.n = n;
        this.bearings = []
        for (var i = 0; i < 2*n; i++) {
            this.bearings.push(i);
        }
    }

    toString() {
        return "Bearings(" + this.bearings.join(",") + ")";
    }

    standardise_index(index) {
        return (index + this.bearings.length) % this.bearings.length;
    }

    get_item(index) {
        var index = this.standardise_index(index);
        return this.bearings[index];
    }

    index_of(item) {
        return this.bearings.indexOf(item)
    }

    collapse_merged() {
        while (true) {
            for (let index = 0; index < this.bearings.length; index++) {
                if (this.get_item(index) % this.n == this.get_item(index+1) % this.n) {
                    if (index + 1 == this.bearings.length) {
                        this.bearings.splice(this.bearings.length -1, 1);
                        this.bearings.splice(0, 1);
                        this.collapse_merged();
                    } else {
                        this.bearings.splice(index, 2);
                        this.collapse_merged();
                    }
                    break;
                }
            }
            break;
        }
    }

    apply_placement(placement) {
        var index_1 = this.index_of(placement.bearing_1)
        var index_2 = this.index_of(placement.bearing_2)
        if (Math.abs(index_1 - index_2) != 1 && Math.abs(index_1 - index_2) != this.bearings.length - 1) {
            console.log('invalid tiling: ' + placement + 'cannot be placed');
            throw EvalError
        }
        [this.bearings[index_1], this.bearings[index_2]] = [this.bearings[index_2], this.bearings[index_1]];

        this.collapse_merged()
        // console.log(this.toString())
    }

    get_valid_placement(already_placed) {
        while (this.bearings.length > 0) {
            var index = Math.floor(Math.random() * this.bearings.length);
            var bearing_1 = this.get_item(index);
            var bearing_2 = this.get_item(index + 1);
            var seen_before = false;
            already_placed.forEach(placement => {
                if (placement.bearing_1 % this.n == bearing_1 % this.n && placement.bearing_2 % this.n == bearing_2 % this.n) {
                    seen_before = true
                }
                if (placement.bearing_1 % this.n == bearing_2 % this.n && placement.bearing_2 % this.n == bearing_1 % this.n) {
                    seen_before = true
                }
            });
            if (seen_before) {
                // console.log("seen before: " + new Placement(this.n, bearing_1, bearing_2).toString())
                continue;
            }
            // console.log("adding: " + new Placement(this.n, bearing_1, bearing_2).toString())
            return new Placement(this.n, bearing_1, bearing_2);
        }
        return null;
    }
}

function get_parallelogram_point(shared_point, other_point_a, other_point_b) {
    var x = other_point_a.x + other_point_b.x - shared_point.x;
    var y = other_point_a.y + other_point_b.y - shared_point.y;
    
    return new Point(x, y)
}

function get_points_from_lines(line_1, line_2) {
    if (Math.abs(line_1.start.x - line_2.start.x) < 2 ** -32 && Math.abs(line_1.start.y - line_2.start.y) < 2 ** -32) {
        return [line_1.end, line_1.start, line_2.end, get_parallelogram_point(line_1.start, line_1.end, line_2.end)]
    }
    if (Math.abs(line_1.start.x - line_2.end.x) < 2 ** -32 && Math.abs(line_1.start.y - line_2.end.y) < 2 ** -32) {
        return [line_1.end, line_1.start, line_2.start, get_parallelogram_point(line_1.start, line_1.end, line_2.start)]
    }
    if (Math.abs(line_1.end.x - line_2.start.x) < 2 ** -32 && Math.abs(line_1.end.y - line_2.start.y) < 2 ** -32) {
        return [line_1.start, line_1.end, line_2.end, get_parallelogram_point(line_1.end, line_1.start, line_2.end)]
    }
    if (Math.abs(line_1.end.x - line_2.end.x) < 2 ** -32 && Math.abs(line_1.end.y - line_2.end.y) < 2 ** -32) {
        return [line_1.start, line_1.end, line_2.start, get_parallelogram_point(line_1.end, line_1.start, line_2.start)]
    }
    console.log('tiling broken: ' + line_1 + ' and ' + line_2 + 'do not share a point');
    throw EvalError
}

function get_perpendicular(point_1, point_2, bearing, n) {
    var difference = new Point(point_2.x - point_1.x, point_2.y - point_1.y);
    var angle = Math.PI * (((bearing/n) + 1) % 2);
    var normal = new Point(Math.cos(angle), Math.sin(angle));
    var dot = difference.x * normal.x + difference.y * normal.y;
    return new Point(dot * normal.x, dot * normal.y);
}

class Tiling {
    constructor(n, face_placements) {
        this.n = n;
        this.placements_to_do = face_placements;
        this.tile_faces = [];
        this.points_array = [[]];
        this.lines_array = [[]];

        // dictionary of bearing to following vertex
        this.bearing_vertices = {};
        for (let index = 0; index < 2 * this.n; index++) {
            this.vertices[index] = new Vertex(n, index);
        }

        // dictionary of bearing to it's edge
        this.bearing_edges = [];
        for (let index = 0; index < 2 * this.n; index++) {
            this.edges[index] = new Edge(n, index);
        }

        // set of tiles which haven't been placed yet
        this.unplaced_tiles = []
        for (var i = 0; i < n; i++) {
            for (var j = i+1; j < n; j++) {
                this.unplaced_tiles.push(new Tile(this.n, i, j));
            }
        }

        // dictionary from bearing to it's open line
        this.bearing_open_lines = [{}]
        for (var bearing in this.edges) {
            var edge = this.edges[bearing];
            this.bearing_open_lines.slice(-1)[0][bearing] = new Line(this.n, edge.start, edge.end, bearing % this.n, bearing >= this.n ? this.n - 1 : 0, true);
        };

        // which bearings are still available and adjacent
        this.bearing_cycle = [new BearingCycle(n)]

        // fill the placements_to_do into their correct places
        while (this.placements_to_do.length > 0) {
            var placement = this.placements_to_do.shift();
            this.add_placement(placement);
        }
        
    }

    toString() {
        return "Tiling(" + this.n + ", " + this.bearing_cycle.slice(-1)[0].join(",") + ")";
    }

    get blueprint() {
        return {
            n: this.n,
            placements: this.tile_faces.map(f => [f.placement.bearing_1, f.placement.bearing_2])
        }
    }

    get vertices() {
        return this.bearing_vertices;
    }

    get points() {
        return this.points_array.slice(-1)[0];
    }

    get edges() {
        return this.bearing_edges;
    }

    get lines() {
        return this.lines_array.slice(-1)[0];
    }

    get faces() {
        return this.tile_faces;
    }

    get curves() {
        var curves = {}
        for (let index = 0; index < this.bearing_open_lines.length; index++) {
            const bearing_lines = this.bearing_open_lines[index];
            for (var bearing in bearing_lines) {
                if (bearing_lines.hasOwnProperty(bearing)) {
                    var line = bearing_lines[bearing]
                    if (!curves.hasOwnProperty(bearing)) {
                        curves[bearing] = []
                    }
                    var centre_x = (line.start.x + line.end.x)/2
                    var centre_y = (line.start.y + line.end.y)/2
                    var previous_point = curves[bearing].slice(-1)[0]
                    if (previous_point == null || centre_x != previous_point.x || centre_y != previous_point.y) {
                        curves[bearing].push(new Point(centre_x, centre_y))
                    }
                }
            }
        }

        return curves;
    }

    bearing_shapes(x_scale, y_scale, edge_colour_type, edge_colour_base) {
        var bearing_shapes = [];
        for (var bearing in this.curves) {
            bearing = parseInt(bearing);
            if (this.curves.hasOwnProperty(bearing)) {
                var points = this.curves[bearing]
                if (points.length < 2) {
                    continue;
                }
                var shape = "M " + [x_scale(points[0].x), y_scale(points[0].y)].join(",")
                for (let index = 1; index < points.length; index++) {
                    const prev = points[index-1];
                    const next = points[index];
                    const centre_x = (prev.x + next.x) / 2
                    const centre_y = (prev.y + next.y) / 2
                    // shape += " T " + [x_scale(centre_x), y_scale(centre_y)].join(",")
                    // shape += " T " + [x_scale(next.x), y_scale(next.y)].join(",")
                    // shape += " C " + [x_scale(centre_x), y_scale(centre_y)].join(" ") + ", " + [x_scale(centre_x), y_scale(centre_y)].join(" ") + ", " + [x_scale(next.x), y_scale(next.y)].join(" ")
                    // shape += " S " + [x_scale(centre_x), y_scale(centre_y)].join(" ") + ", " + [x_scale(next.x), y_scale(next.y)].join(" ")


                    var prev_target = get_perpendicular(prev, next, bearing, this.n);
                    var next_target = get_perpendicular(next, prev, (bearing + this.n) % (2 * this.n), this.n);
                    prev_target = new Point(prev.x + prev_target.x/2, prev.y + prev_target.y/2)
                    next_target = new Point(next.x + next_target.x/2, next.y + next_target.y/2)
                    shape += " C " + [x_scale(prev_target.x), y_scale(prev_target.y)].join(" ") + ", " + [x_scale(next_target.x), y_scale(next_target.y)].join(" ") + ", " + [x_scale(next.x), y_scale(next.y)].join(" ")
                }

                var angle = Math.PI * (((bearing/this.n) + 1) % 2);
                var normal = new Point(Math.cos(angle), Math.sin(angle));
                var prev_target = new Point(points[0].x + normal.x/100, points[0].y + normal.y/100)
                var collapsed_shape = "M " + [x_scale(points[0].x), y_scale(points[0].y)].join(",") + " L " + [x_scale(prev_target.x), y_scale(prev_target.y)].join(" ")

                bearing_shapes.push({
                    bearing: bearing,
                    points: points,
                    shape: shape,
                    collapsed_shape: collapsed_shape,
                    colour: orientation_colour(bearing % this.n, this.n, edge_colour_type, edge_colour_base),
                    key: bearing
                })
            }
        }
        return bearing_shapes;
    }

    get placement_hash() {
        return md5(this.faces.map(d => d.placement))
    }

    get unplaced_table() {
        var face_table = {};
        for (let index = 0; index < this.unplaced_tiles.length; index++) {
            const tile = this.unplaced_tiles[index];
            var start_1 = new Vertex(this.n, tile.orientation_1 - 1);
            var end_1 = new Vertex(this.n, tile.orientation_1);
            var delta_1 = new Point(end_1.x - start_1.x, end_1.y - start_1.y);
            var start_2 = new Vertex(this.n, tile.orientation_2 - 1);
            var end_2 = new Vertex(this.n, tile.orientation_2);
            var delta_2 = new Point(end_2.x - start_2.x, end_2.y - start_2.y);
            var centre = new Point((delta_1.x + delta_2.x)/2, (delta_1.y + delta_2.y)/2)
            var points = [
                new Point(0 - centre.x, 0 - centre.y),
                new Point(delta_1.x - centre.x, delta_1.y - centre.y),
                new Point(delta_1.x + delta_2.x - centre.x, delta_1.y + delta_2.y - centre.y),
                new Point(delta_2.x - centre.x, delta_2.y - centre.y)
            ]
            var face = new Face(this.n, new Placement(this.n, tile.orientation_1, tile.orientation_2), points, false)
            face_table[[tile.orientation_1, tile.orientation_2]] = face
        }

        return face_table
    }

    get unplaced_table_edges() {
        var face_edges = [];
        for (let index = 0; index < this.unplaced_tiles.length; index++) {
            const tile = this.unplaced_tiles[index];
            var start_1 = new Vertex(this.n, tile.orientation_1 - 1);
            var end_1 = new Vertex(this.n, tile.orientation_1);
            var delta_1 = new Point(end_1.x - start_1.x, end_1.y - start_1.y);
            var start_2 = new Vertex(this.n, tile.orientation_2 - 1);
            var end_2 = new Vertex(this.n, tile.orientation_2);
            var delta_2 = new Point(end_2.x - start_2.x, end_2.y - start_2.y);
            var centre = new Point((delta_1.x + delta_2.x)/2, (delta_1.y + delta_2.y)/2)
            face_edges.push(new Line(this.n, new Point(delta_1.x + delta_2.x - centre.x, delta_1.y + delta_2.y - centre.y), new Point(delta_2.x - centre.x, delta_2.y - centre.y), tile.orientation_1, tile.orientation_2, false))
            face_edges.push(new Line(this.n, new Point(0 - centre.x, 0 - centre.y), new Point(delta_2.x - centre.x, delta_2.y - centre.y), tile.orientation_2, tile.orientation_1, false))
        }

        return face_edges
    }

    add_placement(placement) {
        // console.log(placement)
        var already_placed = true;
        for (let index = 0; index < this.unplaced_tiles.length; index++) {
            const tile = this.unplaced_tiles[index];
            if ((tile.orientation_1 == placement.tile.orientation_1 && tile.orientation_2 == placement.tile.orientation_2) || (tile.orientation_1 == placement.tile.orientation_2 && tile.orientation_2 == placement.tile.orientation_1)) {
                already_placed = false;
            }
        }
        if (already_placed) {
            // console.log('invalid tiling: ' + placement + 'already used');
            // console.log(this.unplaced_tiles)
            // throw EvalError
            return null
        }

        // append new data structures
        this.bearing_open_lines.push(clone(this.bearing_open_lines.slice(-1)[0]));
        this.points_array.push(clone(this.points_array.slice(-1)[0]));
        this.lines_array.push(clone(this.lines_array.slice(-1)[0]));
        this.bearing_cycle.push(clone(this.bearing_cycle.slice(-1)[0]));

        // create face and apply data changes
        try {
            this.bearing_cycle.slice(-1)[0].apply_placement(placement)
        } catch (e) {
            this.bearing_open_lines = this.bearing_open_lines.slice(0, -1);
            this.points_array = this.points_array.slice(0, -1);
            this.lines_array = this.lines_array.slice(0, -1);
            this.bearing_cycle = this.bearing_cycle.slice(0, -1);
            return null
        }
        var line_1 = this.bearing_open_lines.slice(-1)[0][placement.bearing_1]
        var line_2 = this.bearing_open_lines.slice(-1)[0][placement.bearing_2]
        var points = get_points_from_lines(line_1, line_2)
        var face = new Face(this.n, placement, points)
        this.tile_faces.push(face);
        points[3].key = this.placement_hash;
        this.points_array.slice(-1)[0].push(points[3]);

        var bearing_1_offset_x = points[2].x - points[1].x
        var bearing_1_offset_y = points[2].y - points[1].y
        var bearing_2_offset_x = points[0].x - points[1].x
        var bearing_2_offset_y = points[0].y - points[1].y

        var start_1 = new Point(this.bearing_open_lines.slice(-1)[0][placement.bearing_1].start.x + bearing_1_offset_x, this.bearing_open_lines.slice(-1)[0][placement.bearing_1].start.y + bearing_1_offset_y)
        var end_1 = new Point(this.bearing_open_lines.slice(-1)[0][placement.bearing_1].end.x + bearing_1_offset_x, this.bearing_open_lines.slice(-1)[0][placement.bearing_1].end.y + bearing_1_offset_y)
        var line_1 = new Line(this.n, start_1, end_1, placement.bearing_1 % this.n, this.bearing_open_lines.slice(-1)[0][placement.bearing_1].step + (placement.bearing_1 >= this.n ? -1 : 1), true)
        this.bearing_open_lines.slice(-1)[0][placement.bearing_1] = line_1
        var opposing_line = this.bearing_open_lines.slice(-1)[0][(placement.bearing_1 + this.n) % (2*this.n)]
        if ((Math.abs(opposing_line.start.x - line_1.start.x) < 2 ** -32 && Math.abs(opposing_line.start.y - line_1.start.y) < 2 ** -32 && Math.abs(opposing_line.end.x - line_1.end.x) < 2 ** -32 && Math.abs(opposing_line.end.y - line_1.end.y) < 2 ** -32)
            || (Math.abs(opposing_line.start.x - line_1.end.x) < 2 ** -32 && Math.abs(opposing_line.start.y - line_1.end.y) < 2 ** -32 && Math.abs(opposing_line.end.x - line_1.start.x) < 2 ** -32 && Math.abs(opposing_line.end.y - line_1.start.y) < 2 ** -32)) {

        } else {
            this.lines_array.slice(-1)[0].push(line_1)
        }

        var start_2 = new Point(this.bearing_open_lines.slice(-1)[0][placement.bearing_2].start.x + bearing_2_offset_x, this.bearing_open_lines.slice(-1)[0][placement.bearing_2].start.y + bearing_2_offset_y)
        var end_2 = new Point(this.bearing_open_lines.slice(-1)[0][placement.bearing_2].end.x + bearing_2_offset_x, this.bearing_open_lines.slice(-1)[0][placement.bearing_2].end.y + bearing_2_offset_y)
        var line_2 = new Line(this.n, start_2, end_2, placement.bearing_2 % this.n, this.bearing_open_lines.slice(-1)[0][placement.bearing_2].step + (placement.bearing_2 >= this.n ? -1 : 1), true)
        this.bearing_open_lines.slice(-1)[0][placement.bearing_2] = line_2
        var opposing_line = this.bearing_open_lines.slice(-1)[0][(placement.bearing_2 + this.n) % (2*this.n)]
        if ((Math.abs(opposing_line.start.x - line_2.start.x) < 2 ** -32 && Math.abs(opposing_line.start.y - line_2.start.y) < 2 ** -32 && Math.abs(opposing_line.end.x - line_2.end.x) < 2 ** -32 && Math.abs(opposing_line.end.y - line_2.end.y) < 2 ** -32)
            || (Math.abs(opposing_line.start.x - line_2.end.x) < 2 ** -32 && Math.abs(opposing_line.start.y - line_2.end.y) < 2 ** -32 && Math.abs(opposing_line.end.x - line_2.start.x) < 2 ** -32 && Math.abs(opposing_line.end.y - line_2.start.y) < 2 ** -32)) {

        } else {
            this.lines_array.slice(-1)[0].push(line_2)
        }

        for (let index = 0; index < this.unplaced_tiles.length; index++) {
            const tile = this.unplaced_tiles[index];
            if (tile.orientation_1 == placement.tile.orientation_1 && tile.orientation_2 == placement.tile.orientation_2) {
                var found_index = index
            }
        }
        this.unplaced_tiles.splice(found_index, 1)
        return face
    }

    add_random_placement() {
        var already_placed = [];
        this.tile_faces.forEach(face => {
            already_placed.push(face.placement)
        });
        var placement = this.bearing_cycle.slice(-1)[0].get_valid_placement(already_placed)
        if (placement == null) {
            return null
        }
        this.add_placement(placement)
        return this
    }

    remove_last_placement() {
        var lost_face = this.tile_faces.splice(this.tile_faces.length - 1)[0];
        this.unplaced_tiles.push(lost_face.placement.tile)
        this.bearing_open_lines = this.bearing_open_lines.slice(0, -1);
        this.points_array = this.points_array.slice(0, -1);
        this.lines_array = this.lines_array.slice(0, -1);
        this.bearing_cycle = this.bearing_cycle.slice(0, -1);

        return this

        // // fill the placements_to_do into their correct places
        // while (this.placements_to_do.length > 0) {
        //     var placement = this.placements_to_do.shift();
        //     this.add_placement(placement);
        // }
        // var new_placements = []
        // for (let index = 0; index < this.faces.length - 1; index++) {
        //     const face = this.faces[index];
        //     new_placements.push(face.placement)
        // }
        // return new Tiling(this.n, new_placements);
    }
    
    remove_tile(tile) {
        var new_placements = [];
        var bearings_found = new Set();
        for (let index = 0; index < this.faces.length; index++) {
            const face = this.faces[index];
            var face_tile = face.placement.tile;
            if ((face_tile.orientation_1 == tile.orientation_1 && face_tile.orientation_2 == tile.orientation_2) || (face_tile.orientation_1 == tile.orientation_2 && face_tile.orientation_2 == tile.orientation_1)) {
                bearings_found.add(face.placement.bearing_1);
                bearings_found.add(face.placement.bearing_2);
            } else {
                if (bearings_found.has(face.placement.bearing_1) || bearings_found.has(face.placement.bearing_2)) {
                    console.log("skipping placement")
                    console.log(face.placement)
                    bearings_found.add(face.placement.bearing_1);
                    bearings_found.add(face.placement.bearing_2);
                } else {
                    new_placements.push(face.placement);
                }
            }
        }
        console.log("removed single tile")
        console.log(new_placements.map(p => "[" + p.bearing_1 + ',' + p.bearing_2 + "]").join(", "))
        return new Tiling(this.n, new_placements);
    }

    decrease_n(new_n) {
        var difference = this.n - new_n;
        var remove_orientations = new Set()
        for (let remove_n = this.n - 1; remove_n > this.n - 1 - difference; remove_n--) {
            remove_orientations.add(remove_n);
        }
        var new_placements = []
        this.tile_faces.forEach(face => {
            var placement = face.placement;
            if (remove_orientations.has(placement.bearing_1 % this.n) || remove_orientations.has(placement.bearing_2 % this.n)) {
            } else {
                var bearing_1 = placement.bearing_1 >= this.n ? placement.bearing_1 - difference : placement.bearing_1;
                var bearing_2 = placement.bearing_2 >= this.n ? placement.bearing_2 - difference : placement.bearing_2;
                var new_placement = new Placement(new_n, bearing_1, bearing_2);
                new_placements.push(new_placement)
            }
        });
        
        return new Tiling(new_n, new_placements);
    }

    increase_n(new_n) {
        var difference = new_n - this.n;
        var new_placements = []
        for (let step_n = new_n; step_n >= this.n + 1; step_n--) {
            for (let index = step_n; index < (2*step_n) - 1; index++) {
                var new_placement = new Placement(new_n, step_n-1, index % (2*step_n) + new_n - step_n)
                new_placements.push(new_placement);
            }
        }
        this.tile_faces.forEach(face => {
            var placement = face.placement;
            var bearing_1 = placement.bearing_1 >= this.n ? placement.bearing_1 + difference : placement.bearing_1;
            var bearing_2 = placement.bearing_2 >= this.n ? placement.bearing_2 + difference : placement.bearing_2;
            var new_placement = new Placement(new_n, bearing_1, bearing_2);
            new_placements.push(new_placement)
        });

        return new Tiling(new_n, new_placements);
    }

    update_n(n) {
        n = parseInt(n);
        if (this.n < n) {
            return this.increase_n(n);
        } else {
            return this.decrease_n(n);
        }
    }
}

function star_tiling(n) {
    placements = [];
    for (let ring = 0; ring < Math.floor(n/2); ring++) {
        for (let i = 0; i < n; i++) {
            placements.push(new Placement(n, 2*i, (1+2*i+2*ring)%(2*n)));
        }
    }

    return placements;
}

// tiling = new Tiling(7,
//     [
//         new Placement(7, 0, 1),
//         new Placement(7, 2, 3),
//         new Placement(7, 5, 6),
//         new Placement(7, 0, 3),
//         new Placement(7, 4, 6),
//         new Placement(7, 2, 6),
//         new Placement(7, 10, 11),
//         // new Placement(7, 10, 12),
//         // new Placement(7, 13, 1),
//         // new Placement(7, 3, 0)
//     ])

// console.log(tiling.faces)
// console.log('done')
