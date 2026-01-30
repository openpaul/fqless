#!/bin/bash
# Generate coverage badges from cargo tarpaulin output

set -e

# Use real grep, not ripgrep alias
GREP=$(command -v /usr/bin/grep || command -v grep)

# Run tarpaulin and capture output
OUTPUT=$(cargo tarpaulin --output-dir target/tarpaulin 2>&1)

# Extract overall coverage percentage
COVERAGE=$(echo "$OUTPUT" | $GREP -oP '^\K[0-9]+\.[0-9]+(?=% coverage)' | head -1)

if [ -z "$COVERAGE" ]; then
    echo "Could not extract coverage percentage"
    exit 1
fi

# Extract per-file coverage
VIEWER_LINES=$(echo "$OUTPUT" | $GREP "src/viewer.rs" | $GREP -oP '\d+/\d+' | head -1)
VIEWER_COVERED=$(echo "$VIEWER_LINES" | cut -d'/' -f1)
VIEWER_TOTAL=$(echo "$VIEWER_LINES" | cut -d'/' -f2)

# Calculate total lines excluding viewer.rs
TOTAL_COVERED=$(echo "$OUTPUT" | $GREP "lines covered" | $GREP -oP '\d+(?=/)')
TOTAL_LINES=$(echo "$OUTPUT" | $GREP "lines covered" | $GREP -oP '(?<=\/)\d+(?= lines)')

# Calculate coverage excluding viewer.rs
if [ -n "$VIEWER_COVERED" ] && [ -n "$VIEWER_TOTAL" ]; then
    NO_TUI_COVERED=$((TOTAL_COVERED - VIEWER_COVERED))
    NO_TUI_TOTAL=$((TOTAL_LINES - VIEWER_TOTAL))
    if [ "$NO_TUI_TOTAL" -gt 0 ]; then
        NO_TUI_COVERAGE=$(awk "BEGIN {printf \"%.2f\", ($NO_TUI_COVERED / $NO_TUI_TOTAL) * 100}")
    else
        NO_TUI_COVERAGE="0.00"
    fi
else
    NO_TUI_COVERAGE="$COVERAGE"
fi

# Function to determine color based on coverage
get_color() {
    local cov=$1
    local cov_int=$(printf "%.0f" "$cov")
    if [ "$cov_int" -ge 80 ]; then
        echo "brightgreen"
    elif [ "$cov_int" -ge 60 ]; then
        echo "green"
    elif [ "$cov_int" -ge 40 ]; then
        echo "yellow"
    elif [ "$cov_int" -ge 20 ]; then
        echo "orange"
    else
        echo "red"
    fi
}

# Function to get hex color
get_hex_color() {
    case $1 in
        brightgreen) echo "#4c1" ;;
        green) echo "#97ca00" ;;
        yellow) echo "#dfb317" ;;
        orange) echo "#fe7d37" ;;
        red) echo "#e05d44" ;;
    esac
}

# Round to integer
COVERAGE_INT=$(printf "%.0f" "$COVERAGE")
NO_TUI_COVERAGE_INT=$(printf "%.0f" "$NO_TUI_COVERAGE")

# Determine colors
COLOR=$(get_color "$COVERAGE")
NO_TUI_COLOR=$(get_color "$NO_TUI_COVERAGE")

# Generate overall coverage badge
# Generate overall coverage badge
cat > coverage-badge.svg << EOF
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="96" height="20">
  <linearGradient id="b" x2="0" y2="100%">
    <stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
    <stop offset="1" stop-opacity=".1"/>
  </linearGradient>
  <clipPath id="a">
    <rect width="96" height="20" rx="3" fill="#fff"/>
  </clipPath>
  <g clip-path="url(#a)">
    <path fill="#555" d="M0 0h59v20H0z"/>
    <path fill="$(get_hex_color $COLOR)" d="M59 0h37v20H59z"/>
    <path fill="url(#b)" d="M0 0h96v20H0z"/>
  </g>
  <g fill="#fff" text-anchor="middle" font-family="DejaVu Sans,Verdana,Geneva,sans-serif" font-size="110">
    <text x="305" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="490">coverage</text>
    <text x="305" y="140" transform="scale(.1)" textLength="490">coverage</text>
    <text x="765" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="270">${COVERAGE_INT}%</text>
    <text x="765" y="140" transform="scale(.1)" textLength="270">${COVERAGE_INT}%</text>
  </g>
</svg>
EOF

# Generate coverage badge excluding TUI
cat > coverage-no-tui-badge.svg << EOF
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="126" height="20">
  <linearGradient id="b" x2="0" y2="100%">
    <stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
    <stop offset="1" stop-opacity=".1"/>
  </linearGradient>
  <clipPath id="a">
    <rect width="126" height="20" rx="3" fill="#fff"/>
  </clipPath>
  <g clip-path="url(#a)">
    <path fill="#555" d="M0 0h89v20H0z"/>
    <path fill="$(get_hex_color $NO_TUI_COLOR)" d="M89 0h37v20H89z"/>
    <path fill="url(#b)" d="M0 0h126v20H0z"/>
  </g>
  <g fill="#fff" text-anchor="middle" font-family="DejaVu Sans,Verdana,Geneva,sans-serif" font-size="110">
    <text x="455" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="790">coverage (no TUI)</text>
    <text x="455" y="140" transform="scale(.1)" textLength="790">coverage (no TUI)</text>
    <text x="1065" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="270">${NO_TUI_COVERAGE_INT}%</text>
    <text x="1065" y="140" transform="scale(.1)" textLength="270">${NO_TUI_COVERAGE_INT}%</text>
  </g>
</svg>
EOF

echo "Coverage badges generated:"
echo "  - coverage-badge.svg (overall: ${COVERAGE}%)"
echo "  - coverage-no-tui-badge.svg (excluding TUI: ${NO_TUI_COVERAGE}%)"