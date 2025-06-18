
use std::fs::File;
use std::io::{BufRead, BufReader};
use image::{RgbImage, Rgb};
use plotters::prelude::*;
use std::collections::HashMap;

/// Parses the input and computes unpaired probabilities
pub fn get_uprobs_from_file(path: &str) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let file = BufReader::new(File::open(path)?);

    let mut uprobs: Vec<Vec<f64>> = Vec::new();
    let mut up: Vec<f64> = Vec::new();

    let mut ltime = String::new();
    let mut llen = 0;

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
            continue;
        }

        let stime = fields[1];
        let occu: f64 = fields[2].parse()?;
        let ss = fields[3];

        if stime == ltime && ss.len() == llen {
            for (j, c) in ss.chars().enumerate() {
                if c == '.' {
                    up[j] += occu;
                }
            }
        } else {
            if !up.is_empty() {
                uprobs.push(up);
            }
            up = vec![0.0; ss.len()];
            for (j, c) in ss.chars().enumerate() {
                if c == '.' {
                    up[j] += occu;
                }
            }
        }

        ltime = stime.to_string();
        llen = ss.len();
    }

    if !up.is_empty() {
        uprobs.push(up);
    }

    Ok(uprobs)
}

/// Writes a heatmap image from a 2D matrix, upscaling each point for visibility.
pub fn write_heatmap_image(
    data: &[Vec<f64>],
    filename: &str,
    scale: u32, // how many pixels per data point (scale = 10 → each dot = 10x10 block)
) -> Result<(), Box<dyn std::error::Error>> {
    let height = data.len();
    let width = data.iter().map(|row| row.len()).max().unwrap_or(0);

    let img_width = width as u32 * scale;
    let img_height = height as u32 * scale;

    let mut img = RgbImage::new(img_width, img_height);

    // Find color range
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    for row in data {
        for &val in row {
            if val > max_val { max_val = val; }
            if val < min_val { min_val = val; }
        }
    }

    for (y, row) in data.iter().enumerate() {
        for (x, &val) in row.iter().enumerate() {
            let norm = (val - min_val) / (max_val - min_val + 1e-9);
            let color = Rgb([
                (norm * 255.0) as u8,
                (norm * 255.0) as u8,
                255u8,
            ]);

            let px = x as u32 * scale;
            let py = y as u32 * scale;

            // Fill a scale×scale block
            for dy in 0..scale {
                for dx in 0..scale {
                    img.put_pixel(px + dx, py + dy, color);
                }
            }
        }
    }

    img.save(filename)?;
    Ok(())
}



/// Render heatmap with flipped Y-axis and correct structure length labeling
pub fn draw_heatmap_plotters(
    data: &[Vec<f64>],
    filename: &str,
    cell_size: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let rows = data.len();
    let cols = data.iter().map(|r| r.len()).max().unwrap_or(0);

    let width = (cols as u32 + 10) * cell_size;
    let height = (rows as u32 + 10) * cell_size;

    let root = BitMapBackend::new(filename, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let (min_val, max_val) = data.iter().flatten().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(min, max), &v| (min.min(v), max.max(v)),
    );

    let color_map = |val: f64| {
        let norm = (val - min_val) / (max_val - min_val + 1e-9);
        RGBColor((norm * 255.0) as u8, (norm * 255.0) as u8, 255u8)
    };

    // 🔁 Map structure length to flipped y-index of first occurrence
    let mut ylabels: HashMap<usize, usize> = HashMap::new();
    for (i, row) in data.iter().enumerate() {
        let len = row.len();
        let flipped_y = rows - i - 1;
        ylabels.entry(len).or_insert(flipped_y);
    }

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .set_label_area_size(LabelAreaPosition::Left, 50)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..cols, 0..(rows as i32))?;

    chart.configure_mesh()
        .disable_mesh()
        .x_labels((cols / 10).max(1))
        .x_label_formatter(&|x| if x % 10 == 0 { format!("{x}") } else { "".into() })
        .x_desc("Sequence Position")
        .y_label_formatter(&|y| {
            let y = *y as usize;
            let label = ylabels.iter().find_map(|(&len, &ypos)| {
                if ypos == y {
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
    for (i, row) in data.iter().enumerate() {
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
    let input_file = "Rnf168-cut.drf";
    let output_file = "heatmap.png";
    let pixel_scale = 2; // ← Change this to make image larger (e.g., 4–20)

    let uprobs = get_uprobs_from_file(input_file)?;
    //write_heatmap_image(&uprobs, output_file, pixel_scale)?;
    draw_heatmap_plotters(&uprobs, output_file, pixel_scale)?;

    println!("Saved heatmap to {}", output_file);
    Ok(())
}

