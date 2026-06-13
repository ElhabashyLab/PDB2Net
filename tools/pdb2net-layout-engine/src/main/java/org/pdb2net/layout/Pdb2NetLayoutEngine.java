package org.pdb2net.layout;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public final class Pdb2NetLayoutEngine {
    private static final double COMPONENT_SPACING = 160.0;
    private static final double MIN_COMPONENT_SIZE = 160.0;
    private static final double OUTPUT_SCALE = 2.0;

    private Pdb2NetLayoutEngine() {
    }

    public static void main(String[] args) throws Exception {
        Map<String, String> cli = parseArgs(args);
        String input = cli.get("input");
        String output = cli.get("output");
        if (input == null || output == null) {
            throw new IllegalArgumentException("Usage: java -jar pdb2net-layout-engine.jar --input job.json --output positions.json");
        }

        Object parsed = new JsonParser(Files.readString(Path.of(input), StandardCharsets.UTF_8)).parse();
        if (!(parsed instanceof Map<?, ?>)) {
            throw new IllegalArgumentException("Input JSON root must be an object");
        }

        LayoutJob job = LayoutJob.fromJson(castMap(parsed));
        Map<String, Point> positions = calculate(job);
        Files.writeString(Path.of(output), writePositions(positions), StandardCharsets.UTF_8);
    }

    private static Map<String, String> parseArgs(String[] args) {
        Map<String, String> values = new HashMap<>();
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if ("--input".equals(arg) || "--output".equals(arg)) {
                if (i + 1 >= args.length) {
                    throw new IllegalArgumentException("Missing value for " + arg);
                }
                values.put(arg.substring(2), args[++i]);
            }
        }
        return values;
    }

    private static Map<String, Point> calculate(LayoutJob job) {
        List<String> nodes = job.nodes;
        int n = nodes.size();
        if (n == 0) {
            return new LinkedHashMap<>();
        }
        if (n == 1) {
            LinkedHashMap<String, Point> single = new LinkedHashMap<>();
            single.put(nodes.get(0), new Point(0.0, 0.0));
            return single;
        }
        if (n == 2 && hasConnectingEdge(nodes.get(0), nodes.get(1), job.edges)) {
            LinkedHashMap<String, Point> pair = new LinkedHashMap<>();
            pair.put(nodes.get(0), new Point(-80.0, 0.0));
            pair.put(nodes.get(1), new Point(80.0, 0.0));
            return pair;
        }

        Map<String, List<String>> adjacency = buildAdjacency(nodes, job.edges);
        List<List<String>> components = connectedComponents(nodes, adjacency);
        components.sort((left, right) -> {
            int sizeCompare = Integer.compare(right.size(), left.size());
            if (sizeCompare != 0) {
                return sizeCompare;
            }
            return left.get(0).compareTo(right.get(0));
        });

        List<ComponentLayout> layouts = new ArrayList<>();
        int index = 0;
        for (List<String> component : components) {
            List<Edge> componentEdges = filterEdges(component, job.edges);
            long seed = deterministicSeed(job.networkTitle + "|" + index + "|" + component.toString());
            Map<String, Point> local = layoutComponent(component, componentEdges, job.params, seed);
            layouts.add(new ComponentLayout(component, normalize(local)));
            index++;
        }

        return packComponents(layouts, nodes);
    }

    private static boolean hasConnectingEdge(String a, String b, List<Edge> edges) {
        for (Edge edge : edges) {
            if ((edge.source.equals(a) && edge.target.equals(b)) || (edge.source.equals(b) && edge.target.equals(a))) {
                return true;
            }
        }
        return false;
    }

    private static Map<String, List<String>> buildAdjacency(List<String> nodes, List<Edge> edges) {
        Map<String, List<String>> adjacency = new LinkedHashMap<>();
        for (String node : nodes) {
            adjacency.put(node, new ArrayList<>());
        }
        for (Edge edge : edges) {
            if (edge.source.equals(edge.target)) {
                continue;
            }
            if (!adjacency.containsKey(edge.source) || !adjacency.containsKey(edge.target)) {
                continue;
            }
            adjacency.get(edge.source).add(edge.target);
            adjacency.get(edge.target).add(edge.source);
        }
        for (List<String> neighbors : adjacency.values()) {
            Collections.sort(neighbors);
        }
        return adjacency;
    }

    private static List<List<String>> connectedComponents(List<String> nodes, Map<String, List<String>> adjacency) {
        List<List<String>> components = new ArrayList<>();
        Set<String> visited = new HashSet<>();
        for (String start : nodes) {
            if (visited.contains(start)) {
                continue;
            }
            List<String> component = new ArrayList<>();
            ArrayDeque<String> queue = new ArrayDeque<>();
            visited.add(start);
            queue.add(start);
            while (!queue.isEmpty()) {
                String current = queue.removeFirst();
                component.add(current);
                for (String neighbor : adjacency.getOrDefault(current, List.of())) {
                    if (visited.add(neighbor)) {
                        queue.addLast(neighbor);
                    }
                }
            }
            Collections.sort(component);
            components.add(component);
        }
        return components;
    }

    private static List<Edge> filterEdges(List<String> component, List<Edge> edges) {
        Set<String> members = new HashSet<>(component);
        List<Edge> filtered = new ArrayList<>();
        Set<String> seen = new LinkedHashSet<>();
        for (Edge edge : edges) {
            if (edge.source.equals(edge.target)) {
                continue;
            }
            if (!members.contains(edge.source) || !members.contains(edge.target)) {
                continue;
            }
            String a = edge.source.compareTo(edge.target) <= 0 ? edge.source : edge.target;
            String b = edge.source.compareTo(edge.target) <= 0 ? edge.target : edge.source;
            String key = a + "\u0000" + b;
            if (seen.add(key)) {
                filtered.add(new Edge(a, b));
            }
        }
        return filtered;
    }

    private static Map<String, Point> layoutComponent(List<String> nodes, List<Edge> edges, LayoutParams params, long seed) {
        LinkedHashMap<String, Point> positions = new LinkedHashMap<>();
        int n = nodes.size();
        if (n == 0) {
            return positions;
        }
        if (n == 1) {
            positions.put(nodes.get(0), new Point(0.0, 0.0));
            return positions;
        }
        if (n == 2 && !edges.isEmpty()) {
            positions.put(nodes.get(0), new Point(-80.0, 0.0));
            positions.put(nodes.get(1), new Point(80.0, 0.0));
            return positions;
        }

        Random random = new Random(seed);
        double initRadius = Math.max(80.0, Math.sqrt(n) * params.springLength * 0.8);
        for (int i = 0; i < n; i++) {
            double angle = (2.0 * Math.PI * i / n) + (random.nextDouble() - 0.5) * 0.12;
            double radius = initRadius * (0.75 + random.nextDouble() * 0.25);
            positions.put(nodes.get(i), new Point(Math.cos(angle) * radius, Math.sin(angle) * radius));
        }

        double velocityDamping = 0.78;
        double timeStep = 0.35;
        Map<String, Point> velocity = new HashMap<>();
        for (String node : nodes) {
            velocity.put(node, new Point(0.0, 0.0));
        }

        for (int iteration = 0; iteration < params.iterations; iteration++) {
            Map<String, Point> forces = new HashMap<>();
            for (String node : nodes) {
                forces.put(node, new Point(0.0, 0.0));
            }

            for (int i = 0; i < n; i++) {
                String a = nodes.get(i);
                Point pa = positions.get(a);
                for (int j = i + 1; j < n; j++) {
                    String b = nodes.get(j);
                    Point pb = positions.get(b);
                    double dx = pa.x - pb.x;
                    double dy = pa.y - pb.y;
                    double distanceSquared = Math.max(25.0, dx * dx + dy * dy);
                    double distance = Math.sqrt(distanceSquared);
                    double force = (params.nodeMass * params.nodeMass * 850.0) / distanceSquared;
                    double fx = (dx / distance) * force;
                    double fy = (dy / distance) * force;
                    forces.get(a).add(fx, fy);
                    forces.get(b).add(-fx, -fy);
                }
            }

            for (Edge edge : edges) {
                Point source = positions.get(edge.source);
                Point target = positions.get(edge.target);
                if (source == null || target == null) {
                    continue;
                }
                double dx = target.x - source.x;
                double dy = target.y - source.y;
                double distance = Math.max(0.001, Math.sqrt(dx * dx + dy * dy));
                double displacement = distance - params.springLength;
                double force = params.springCoefficient * displacement * 10000.0;
                double fx = (dx / distance) * force;
                double fy = (dy / distance) * force;
                forces.get(edge.source).add(fx, fy);
                forces.get(edge.target).add(-fx, -fy);
            }

            for (String node : nodes) {
                Point v = velocity.get(node);
                Point f = forces.get(node);
                v.x = (v.x + f.x * timeStep) * velocityDamping;
                v.y = (v.y + f.y * timeStep) * velocityDamping;
                double speed = Math.sqrt(v.x * v.x + v.y * v.y);
                double maxSpeed = 35.0;
                if (speed > maxSpeed) {
                    v.x = (v.x / speed) * maxSpeed;
                    v.y = (v.y / speed) * maxSpeed;
                }
                positions.get(node).add(v.x, v.y);
            }
        }

        return positions;
    }

    private static Map<String, Point> normalize(Map<String, Point> positions) {
        if (positions.isEmpty()) {
            return positions;
        }
        double cx = 0.0;
        double cy = 0.0;
        for (Point point : positions.values()) {
            cx += point.x;
            cy += point.y;
        }
        cx /= positions.size();
        cy /= positions.size();
        LinkedHashMap<String, Point> normalized = new LinkedHashMap<>();
        for (Map.Entry<String, Point> entry : positions.entrySet()) {
            normalized.put(entry.getKey(), new Point(entry.getValue().x - cx, entry.getValue().y - cy));
        }
        return normalized;
    }

    private static Map<String, Point> packComponents(List<ComponentLayout> layouts, List<String> originalNodeOrder) {
        List<PackedComponent> packed = new ArrayList<>();
        for (ComponentLayout layout : layouts) {
            Bounds bounds = Bounds.from(layout.positions);
            double width = Math.max(MIN_COMPONENT_SIZE, bounds.width());
            double height = Math.max(MIN_COMPONENT_SIZE, bounds.height());
            packed.add(new PackedComponent(layout, bounds, width, height));
        }

        double totalArea = 0.0;
        for (PackedComponent component : packed) {
            totalArea += component.width * component.height;
        }
        double targetWidth = Math.max(MIN_COMPONENT_SIZE, Math.sqrt(totalArea) * 1.4);

        double cursorX = 0.0;
        double cursorY = 0.0;
        double rowHeight = 0.0;
        for (PackedComponent component : packed) {
            if (cursorX > 0.0 && cursorX + component.width > targetWidth) {
                cursorX = 0.0;
                cursorY += rowHeight + COMPONENT_SPACING;
                rowHeight = 0.0;
            }
            component.offsetX = cursorX - component.bounds.minX;
            component.offsetY = cursorY - component.bounds.minY;
            cursorX += component.width + COMPONENT_SPACING;
            rowHeight = Math.max(rowHeight, component.height);
        }

        LinkedHashMap<String, Point> finalPositions = new LinkedHashMap<>();
        for (PackedComponent component : packed) {
            for (Map.Entry<String, Point> entry : component.layout.positions.entrySet()) {
                finalPositions.put(
                    entry.getKey(),
                    new Point(entry.getValue().x + component.offsetX, entry.getValue().y + component.offsetY)
                );
            }
        }

        Bounds finalBounds = Bounds.from(finalPositions);
        double centerX = finalBounds.minX + finalBounds.width() / 2.0;
        double centerY = finalBounds.minY + finalBounds.height() / 2.0;
        LinkedHashMap<String, Point> ordered = new LinkedHashMap<>();
        for (String node : originalNodeOrder) {
            Point point = finalPositions.getOrDefault(node, new Point(0.0, 0.0));
            ordered.put(node, new Point((point.x - centerX) * OUTPUT_SCALE, (point.y - centerY) * OUTPUT_SCALE));
        }
        return ordered;
    }

    private static long deterministicSeed(String text) {
        long hash = 1125899906842597L;
        for (int i = 0; i < text.length(); i++) {
            hash = 31 * hash + text.charAt(i);
        }
        return hash;
    }

    private static String writePositions(Map<String, Point> positions) {
        StringBuilder builder = new StringBuilder();
        builder.append("{\n  \"positions\": {");
        int index = 0;
        for (Map.Entry<String, Point> entry : positions.entrySet()) {
            if (index > 0) {
                builder.append(",");
            }
            builder.append("\n    ");
            appendJsonString(builder, entry.getKey());
            builder.append(": {\"x\": ");
            builder.append(formatDouble(entry.getValue().x));
            builder.append(", \"y\": ");
            builder.append(formatDouble(entry.getValue().y));
            builder.append("}");
            index++;
        }
        if (!positions.isEmpty()) {
            builder.append("\n  ");
        }
        builder.append("}\n}\n");
        return builder.toString();
    }

    private static String formatDouble(double value) {
        if (!Double.isFinite(value)) {
            return "0.0";
        }
        return String.format(Locale.ROOT, "%.6f", value);
    }

    private static void appendJsonString(StringBuilder builder, String value) {
        builder.append('"');
        for (int i = 0; i < value.length(); i++) {
            char ch = value.charAt(i);
            switch (ch) {
                case '\\':
                    builder.append("\\\\");
                    break;
                case '"':
                    builder.append("\\\"");
                    break;
                case '\n':
                    builder.append("\\n");
                    break;
                case '\r':
                    builder.append("\\r");
                    break;
                case '\t':
                    builder.append("\\t");
                    break;
                default:
                    if (ch < 0x20) {
                        builder.append(String.format(Locale.ROOT, "\\u%04x", (int) ch));
                    } else {
                        builder.append(ch);
                    }
            }
        }
        builder.append('"');
    }

    @SuppressWarnings("unchecked")
    private static Map<String, Object> castMap(Object value) {
        return (Map<String, Object>) value;
    }

    @SuppressWarnings("unchecked")
    private static List<Object> castList(Object value) {
        return (List<Object>) value;
    }

    private static String stringValue(Object value, String fallback) {
        return value == null ? fallback : String.valueOf(value);
    }

    private static double doubleValue(Map<String, Object> map, String key, double fallback) {
        Object value = map.get(key);
        if (value instanceof Number number) {
            return number.doubleValue();
        }
        if (value != null) {
            try {
                return Double.parseDouble(String.valueOf(value));
            } catch (NumberFormatException ignored) {
                return fallback;
            }
        }
        return fallback;
    }

    private static int intValue(Map<String, Object> map, String key, int fallback) {
        Object value = map.get(key);
        if (value instanceof Number number) {
            return number.intValue();
        }
        if (value != null) {
            try {
                return Integer.parseInt(String.valueOf(value));
            } catch (NumberFormatException ignored) {
                return fallback;
            }
        }
        return fallback;
    }

    private record Edge(String source, String target) {
    }

    private static final class LayoutJob {
        final String networkTitle;
        final List<String> nodes;
        final List<Edge> edges;
        final LayoutParams params;

        LayoutJob(String networkTitle, List<String> nodes, List<Edge> edges, LayoutParams params) {
            this.networkTitle = networkTitle;
            this.nodes = nodes;
            this.edges = edges;
            this.params = params;
        }

        static LayoutJob fromJson(Map<String, Object> root) {
            String title = stringValue(root.get("network_title"), "");
            List<String> nodes = new ArrayList<>();
            Object rawNodes = root.get("nodes");
            if (rawNodes instanceof List<?>) {
                for (Object item : castList(rawNodes)) {
                    if (item instanceof Map<?, ?> nodeMap) {
                        String id = stringValue(castMap(nodeMap).get("id"), "");
                        if (!id.isBlank()) {
                            nodes.add(id);
                        }
                    }
                }
            }

            LinkedHashSet<String> nodeSet = new LinkedHashSet<>(nodes);
            List<Edge> edges = new ArrayList<>();
            Object rawEdges = root.get("edges");
            if (rawEdges instanceof List<?>) {
                for (Object item : castList(rawEdges)) {
                    if (item instanceof Map<?, ?> edgeMap) {
                        Map<String, Object> map = castMap(edgeMap);
                        String source = stringValue(map.get("source"), "");
                        String target = stringValue(map.get("target"), "");
                        if (!source.isBlank() && !target.isBlank() && nodeSet.contains(source) && nodeSet.contains(target)) {
                            edges.add(new Edge(source, target));
                        }
                    }
                }
            }

            LayoutParams params = LayoutParams.defaults();
            Object rawLayout = root.get("layout");
            if (rawLayout instanceof Map<?, ?> layoutMap) {
                Map<String, Object> layout = castMap(layoutMap);
                params = new LayoutParams(
                    intValue(layout, "iterations", params.iterations),
                    doubleValue(layout, "spring_coefficient", params.springCoefficient),
                    doubleValue(layout, "spring_length", params.springLength),
                    doubleValue(layout, "node_mass", params.nodeMass)
                );
            }
            return new LayoutJob(title, nodes, edges, params);
        }
    }

    private static final class LayoutParams {
        final int iterations;
        final double springCoefficient;
        final double springLength;
        final double nodeMass;

        LayoutParams(int iterations, double springCoefficient, double springLength, double nodeMass) {
            this.iterations = Math.max(1, iterations);
            this.springCoefficient = springCoefficient;
            this.springLength = Math.max(1.0, springLength);
            this.nodeMass = Math.max(0.1, nodeMass);
        }

        static LayoutParams defaults() {
            return new LayoutParams(100, 0.0001, 50.0, 3.0);
        }
    }

    private static final class Point {
        double x;
        double y;

        Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        void add(double dx, double dy) {
            this.x += dx;
            this.y += dy;
        }
    }

    private static final class Bounds {
        final double minX;
        final double minY;
        final double maxX;
        final double maxY;

        Bounds(double minX, double minY, double maxX, double maxY) {
            this.minX = minX;
            this.minY = minY;
            this.maxX = maxX;
            this.maxY = maxY;
        }

        double width() {
            return maxX - minX;
        }

        double height() {
            return maxY - minY;
        }

        static Bounds from(Map<String, Point> positions) {
            if (positions.isEmpty()) {
                return new Bounds(0.0, 0.0, 0.0, 0.0);
            }
            double minX = Double.POSITIVE_INFINITY;
            double minY = Double.POSITIVE_INFINITY;
            double maxX = Double.NEGATIVE_INFINITY;
            double maxY = Double.NEGATIVE_INFINITY;
            for (Point point : positions.values()) {
                minX = Math.min(minX, point.x);
                minY = Math.min(minY, point.y);
                maxX = Math.max(maxX, point.x);
                maxY = Math.max(maxY, point.y);
            }
            return new Bounds(minX, minY, maxX, maxY);
        }
    }

    private static final class ComponentLayout {
        final List<String> nodes;
        final Map<String, Point> positions;

        ComponentLayout(List<String> nodes, Map<String, Point> positions) {
            this.nodes = nodes;
            this.positions = positions;
        }
    }

    private static final class PackedComponent {
        final ComponentLayout layout;
        final Bounds bounds;
        final double width;
        final double height;
        double offsetX;
        double offsetY;

        PackedComponent(ComponentLayout layout, Bounds bounds, double width, double height) {
            this.layout = layout;
            this.bounds = bounds;
            this.width = width;
            this.height = height;
        }
    }

    private static final class JsonParser {
        private final String text;
        private int index;

        JsonParser(String text) {
            this.text = text;
        }

        Object parse() {
            Object value = parseValue();
            skipWhitespace();
            if (index != text.length()) {
                throw new IllegalArgumentException("Unexpected trailing JSON at index " + index);
            }
            return value;
        }

        private Object parseValue() {
            skipWhitespace();
            if (index >= text.length()) {
                throw new IllegalArgumentException("Unexpected end of JSON");
            }
            char ch = text.charAt(index);
            if (ch == '{') {
                return parseObject();
            }
            if (ch == '[') {
                return parseArray();
            }
            if (ch == '"') {
                return parseString();
            }
            if (ch == 't' && text.startsWith("true", index)) {
                index += 4;
                return Boolean.TRUE;
            }
            if (ch == 'f' && text.startsWith("false", index)) {
                index += 5;
                return Boolean.FALSE;
            }
            if (ch == 'n' && text.startsWith("null", index)) {
                index += 4;
                return null;
            }
            return parseNumber();
        }

        private Map<String, Object> parseObject() {
            expect('{');
            LinkedHashMap<String, Object> map = new LinkedHashMap<>();
            skipWhitespace();
            if (peek('}')) {
                index++;
                return map;
            }
            while (true) {
                skipWhitespace();
                String key = parseString();
                skipWhitespace();
                expect(':');
                Object value = parseValue();
                map.put(key, value);
                skipWhitespace();
                if (peek('}')) {
                    index++;
                    return map;
                }
                expect(',');
            }
        }

        private List<Object> parseArray() {
            expect('[');
            List<Object> list = new ArrayList<>();
            skipWhitespace();
            if (peek(']')) {
                index++;
                return list;
            }
            while (true) {
                list.add(parseValue());
                skipWhitespace();
                if (peek(']')) {
                    index++;
                    return list;
                }
                expect(',');
            }
        }

        private String parseString() {
            expect('"');
            StringBuilder builder = new StringBuilder();
            while (index < text.length()) {
                char ch = text.charAt(index++);
                if (ch == '"') {
                    return builder.toString();
                }
                if (ch == '\\') {
                    if (index >= text.length()) {
                        throw new IllegalArgumentException("Unterminated JSON escape");
                    }
                    char esc = text.charAt(index++);
                    switch (esc) {
                        case '"':
                        case '\\':
                        case '/':
                            builder.append(esc);
                            break;
                        case 'b':
                            builder.append('\b');
                            break;
                        case 'f':
                            builder.append('\f');
                            break;
                        case 'n':
                            builder.append('\n');
                            break;
                        case 'r':
                            builder.append('\r');
                            break;
                        case 't':
                            builder.append('\t');
                            break;
                        case 'u':
                            if (index + 4 > text.length()) {
                                throw new IllegalArgumentException("Invalid unicode escape");
                            }
                            String hex = text.substring(index, index + 4);
                            builder.append((char) Integer.parseInt(hex, 16));
                            index += 4;
                            break;
                        default:
                            throw new IllegalArgumentException("Invalid JSON escape: " + esc);
                    }
                } else {
                    builder.append(ch);
                }
            }
            throw new IllegalArgumentException("Unterminated JSON string");
        }

        private Number parseNumber() {
            int start = index;
            if (peek('-')) {
                index++;
            }
            while (index < text.length() && Character.isDigit(text.charAt(index))) {
                index++;
            }
            if (peek('.')) {
                index++;
                while (index < text.length() && Character.isDigit(text.charAt(index))) {
                    index++;
                }
            }
            if (index < text.length() && (text.charAt(index) == 'e' || text.charAt(index) == 'E')) {
                index++;
                if (index < text.length() && (text.charAt(index) == '+' || text.charAt(index) == '-')) {
                    index++;
                }
                while (index < text.length() && Character.isDigit(text.charAt(index))) {
                    index++;
                }
            }
            String raw = text.substring(start, index);
            if (raw.isBlank()) {
                throw new IllegalArgumentException("Expected JSON number at index " + start);
            }
            if (raw.contains(".") || raw.contains("e") || raw.contains("E")) {
                return Double.parseDouble(raw);
            }
            try {
                return Integer.parseInt(raw);
            } catch (NumberFormatException exc) {
                return Long.parseLong(raw);
            }
        }

        private void expect(char expected) {
            skipWhitespace();
            if (index >= text.length() || text.charAt(index) != expected) {
                throw new IllegalArgumentException("Expected '" + expected + "' at index " + index);
            }
            index++;
        }

        private boolean peek(char expected) {
            return index < text.length() && text.charAt(index) == expected;
        }

        private void skipWhitespace() {
            while (index < text.length()) {
                char ch = text.charAt(index);
                if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') {
                    index++;
                } else {
                    return;
                }
            }
        }
    }
}
