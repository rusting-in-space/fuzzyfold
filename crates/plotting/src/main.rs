
use std::fs::File;
use std::io::{BufRead, BufReader};
use plotters::prelude::*;
use std::collections::HashMap;

/// Parses the input and computes unpaired probabilities
pub fn get_uprobs_from_drf_file(path: &str) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let file = BufReader::new(File::open(path)?);

    let mut uprobs: Vec<Vec<f64>> = Vec::new();
    let mut up: Vec<f64> = Vec::new();

    let mut ltime = String::new();
    let mut llen = 0;

    fn up_add_occu(up: &mut Vec<f64>, dotbr: &str, occu: f64) {
        for (j, c) in dotbr.chars().enumerate() {
            if c == '.' {
                up[j] += occu;
            }
        }
    }

    for (i, line) in file.lines().enumerate() {
        let line = line?;
        if i == 0 {
            if line.trim() != "id time occupancy structure energy" {
                return Err("Header does not match".into());
            }
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 {
            //TODO?
            continue;
        }

        let stime = fields[1];
        let occu: f64 = fields[2].parse()?;
        let dotbr = fields[3];

        if stime == ltime && dotbr.len() == llen {
            up_add_occu(&mut up, &dotbr, occu);
        } else {
            if !up.is_empty() {
                uprobs.push(up);
            }
            up = vec![0.0; dotbr.len()];
            up_add_occu(&mut up, &dotbr, occu);
        }

        ltime = stime.to_string();
        llen = dotbr.len();
    }

    if !up.is_empty() {
        uprobs.push(up);
    }

    Ok(uprobs)
}

/// Render heatmap with flipped Y-axis and correct structure length labeling
pub fn draw_accessibility_heatmap(
    uprobs: &[Vec<f64>],
    filename: &str,
    cell_size: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let rows = uprobs.len();
    let cols = uprobs.iter().map(|r| r.len()).max().unwrap_or(0);

    let width = (cols as u32 + 10) * cell_size;
    let height = (rows as u32 + 10) * cell_size;

    let root = BitMapBackend::new(filename, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let (min_val, max_val) = uprobs.iter().flatten().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(min, max), &v| (min.min(v), max.max(v)),
    );

    let color_map = |val: f64| {
        let norm = (val - min_val) / (max_val - min_val + 1e-9);
        RGBColor((norm * 255.0) as u8, (norm * 255.0) as u8, 255u8)
    };

    // 🔁 Map structure length to flipped y-index of first occurrence
    let mut ylabels: HashMap<usize, usize> = HashMap::new();
    for (i, row) in uprobs.iter().enumerate() {
        let len = row.len();
        let flipped_y = rows - i - 1;
        ylabels.entry(len).or_insert(flipped_y);
    }

    let mut chart = ChartBuilder::on(&root)
        .caption("Simulated accessibilities", ("sans-serif", (5).percent_width()))
        .margin(20)
        .set_label_area_size(LabelAreaPosition::Left, 50)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..cols, 0..(rows as i32))?;

    chart.configure_mesh()
        .disable_mesh()
        .x_labels((cols / 10).max(1))
        .x_label_formatter(&|x| if x % 10 == 0 { format!("{x}") } else { "".into() })
        .x_desc("Sequence Position")
        .y_labels(rows)
        .y_label_formatter(&|y| {
            let y = *y as usize;
            let label = ylabels.iter().find_map(|(&len, &ypos)| {
                if ypos == y && len % 10 == 0 {
                    Some(format!("{len}"))
                } else {
                    None
                }
            });
            label.unwrap_or_else(|| "".into())
        })
        .y_desc("Structure Length")
        .draw()?;

    // draw rectangles with flipped y-axis
    for (i, row) in uprobs.iter().enumerate() {
        let y = (rows - i - 1) as i32;
        for (x, &val) in row.iter().enumerate() {
            let color = color_map(val).filled();
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x, y), (x + 1, y + 1)],
                color,
            )))?;
        }
    }

    Ok(())
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input_file = "Rnf.drf";
    let output_file = "heatmap.png";
    let pixel_scale = 1; // ← Change this to make image larger (e.g., 4–20)

    let uprobs = get_uprobs_from_drf_file(input_file)?;
    draw_accessibility_heatmap(&uprobs, output_file, pixel_scale)?;

    println!("Saved heatmap to {}", output_file);
    Ok(())
}

