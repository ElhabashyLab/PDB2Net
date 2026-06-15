#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
CLASS_DIR="$BUILD_DIR/classes"
JAR_PATH="$BUILD_DIR/pdb2net-layout-engine.jar"

rm -rf "$CLASS_DIR"
mkdir -p "$CLASS_DIR"

javac -encoding UTF-8 -d "$CLASS_DIR" "$ROOT_DIR/src/main/java/org/pdb2net/layout/Pdb2NetLayoutEngine.java"
jar --create --file "$JAR_PATH" --main-class org.pdb2net.layout.Pdb2NetLayoutEngine -C "$CLASS_DIR" .

echo "$JAR_PATH"
