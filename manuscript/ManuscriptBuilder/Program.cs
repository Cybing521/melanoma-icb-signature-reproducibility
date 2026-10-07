// MDPI 稿件生成器（后端 S：直接 OpenXML）
//
// 为什么不用 CLI：内置 create 命令的 content-json 只支持 heading / paragraph /
// pagebreak，装不下论文必需的表格、图与行号。
//
// 为什么把它留在仓库里：本项目的核心主张是「一切可追溯」。
// 稿件正文里每一个数字都必须能追回某个脚本的某次运行，
// 因此稿件的生成过程本身也必须是可复现的，而不是一次性的手工排版。
//
// 体例依据 MDPI 投稿要求：
//   - A4，双倍行距，Times New Roman 12 pt
//   - 连续行号（MDPI 审稿硬性要求）
//   - 章节编号 1. / 2.1. / 2.1.1.
//   - 表格题注在表上方、图题注在图下方
//
// 用法：dotnet run --project manuscript/ManuscriptBuilder -- <content.json> <out.docx> <baseDir>

using System.Text.Json;
using DocumentFormat.OpenXml;
using DocumentFormat.OpenXml.Packaging;
using DocumentFormat.OpenXml.Wordprocessing;

namespace ManuscriptBuilder;

internal static class Program
{
    // MDPI 投稿稿：Times New Roman 12 pt
    private const string FontName = "Times New Roman";
    private const int BodyHalfPt = 24;          // 12 pt
    private const int LineTwips = 480;          // 双倍行距
    private const int DxaPerMm = 56;            // 1 mm = 56 DXA

    private static int Main(string[] args)
    {
        if (args.Length < 3)
        {
            Console.Error.WriteLine("usage: ManuscriptBuilder <content.json> <out.docx> <baseDir>");
            return 2;
        }
        var contentPath = args[0];
        var outPath = args[1];
        var baseDir = args[2];

        using var doc = JsonDocument.Parse(File.ReadAllText(contentPath));
        var root = doc.RootElement;

        if (File.Exists(outPath)) File.Delete(outPath);
        using var w = WordprocessingDocument.Create(outPath, WordprocessingDocumentType.Document);
        var main = w.AddMainDocumentPart();
        main.Document = new Document();
        var body = main.Document.AppendChild(new Body());

        AddStyles(main);

        if (root.TryGetProperty("meta", out var meta))
        {
            AddFrontMatter(body, meta);
        }

        foreach (var blk in root.GetProperty("blocks").EnumerateArray())
        {
            switch (blk.GetProperty("type").GetString())
            {
                case "heading": AddHeading(body, blk); break;
                case "paragraph": AddParagraph(body, blk); break;
                case "bullet": AddBullet(body, blk); break;
                case "table": AddTable(body, blk, baseDir); break;
                case "figure": AddFigure(body, blk, baseDir); break;
                case "pagebreak": body.AppendChild(new Paragraph(new Run(new Break { Type = BreakValues.Page }))); break;
                default: Console.Error.WriteLine($"warn: unknown block type {blk.GetProperty("type")}"); break;
            }
        }

        body.AppendChild(BuildSectionProperties());
        main.Document.Save();

        // 合并同格式相邻 run，减小体积
        var cli = Environment.GetEnvironmentVariable("DOCX_CLI");
        Console.WriteLine($"written: {outPath}");
        return 0;
    }

    // ── 节属性：A4 + 页边距 + 连续行号 ────────────────────────────
    private static SectionProperties BuildSectionProperties()
    {
        var sect = new SectionProperties();
        sect.Append(new PageSize { Width = 11906U, Height = 16838U });   // A4
        sect.Append(new PageMargin
        {
            Top = 1440,
            Right = 1440U,
            Bottom = 1440,
            Left = 1440U,
            Header = 720U,
            Footer = 720U,
            Gutter = 0U
        });
        // MDPI 审稿要求连续行号
        sect.Append(new LineNumberType
        {
            CountBy = (short)1,
            Restart = LineNumberRestartValues.Continuous,
            Distance = "360"
        });
        return sect;
    }

    // ── 样式表 ────────────────────────────────────────────────────
    private static void AddStyles(MainDocumentPart main)
    {
        var part = main.AddNewPart<StyleDefinitionsPart>();
        var styles = new Styles();

        styles.Append(new DocDefaults(
            new RunPropertiesDefault(new RunPropertiesBaseStyle(
                new RunFonts { Ascii = FontName, HighAnsi = FontName, ComplexScript = FontName },
                new FontSize { Val = BodyHalfPt.ToString() },
                new FontSizeComplexScript { Val = BodyHalfPt.ToString() })),
            new ParagraphPropertiesDefault(new ParagraphPropertiesBaseStyle(
                new SpacingBetweenLines { After = "0", Line = LineTwips.ToString(), LineRule = LineSpacingRuleValues.Auto }))));

        styles.Append(MakeStyle("Normal", "Normal", isDefault: true, sizeHalfPt: BodyHalfPt, bold: false, spaceBefore: 0));

        // 标题：MDPI 用无衬线粗体；此处保持 Times 以符合多数投稿模板
        styles.Append(MakeStyle("Heading1", "heading 1", sizeHalfPt: 28, bold: true, spaceBefore: 240, outline: 0));
        styles.Append(MakeStyle("Heading2", "heading 2", sizeHalfPt: 26, bold: true, spaceBefore: 200, outline: 1));
        styles.Append(MakeStyle("Heading3", "heading 3", sizeHalfPt: 24, bold: true, spaceBefore: 160, outline: 2));
        styles.Append(MakeStyle("Title", "Title", sizeHalfPt: 32, bold: true, spaceBefore: 0));
        styles.Append(MakeStyle("Caption", "caption", sizeHalfPt: 20, bold: false, spaceBefore: 80));
        styles.Append(MakeStyle("TableText", "Table Text", sizeHalfPt: 16, bold: false, spaceBefore: 20));

        part.Styles = styles;
        part.Styles.Save();
    }

    private static Style MakeStyle(string id, string name, bool isDefault = false,
                                  int sizeHalfPt = 24, bool bold = false,
                                  int spaceBefore = 0, int? outline = null)
    {
        var st = new Style
        {
            Type = StyleValues.Paragraph,
            StyleId = id,
            CustomStyle = !isDefault
        };
        if (isDefault) st.SetAttribute(new OpenXmlAttribute("", "default", "", "1"));
        st.Append(new StyleName { Val = name });
        if (id != "Normal") st.Append(new BasedOn { Val = "Normal" });

        var pPr = new StyleParagraphProperties(
            new SpacingBetweenLines { Before = spaceBefore.ToString(), After = "0", Line = LineTwips.ToString(), LineRule = LineSpacingRuleValues.Auto },
            new Justification { Val = JustificationValues.Both });
        if (outline.HasValue) pPr.Append(new OutlineLevel { Val = outline.Value });
        st.Append(pPr);

        var rPr = new StyleRunProperties(
            new RunFonts { Ascii = FontName, HighAnsi = FontName, ComplexScript = FontName },
            new FontSize { Val = sizeHalfPt.ToString() },
            new FontSizeComplexScript { Val = sizeHalfPt.ToString() });
        if (bold) rPr.Append(new Bold());
        st.Append(rPr);
        return st;
    }

    // ── 题名与作者块 ──────────────────────────────────────────────
    private static void AddFrontMatter(Body body, JsonElement meta)
    {
        if (meta.TryGetProperty("title", out var t))
        {
            body.AppendChild(new Paragraph(
                new ParagraphProperties(new ParagraphStyleId { Val = "Title" },
                                        new Justification { Val = JustificationValues.Center }),
                new Run(new Text(t.GetString() ?? ""))));
        }
        foreach (var key in new[] { "authors", "affiliations", "correspondence" })
        {
            if (meta.TryGetProperty(key, out var v))
            {
                body.AppendChild(new Paragraph(
                    new ParagraphProperties(new Justification { Val = JustificationValues.Center },
                                            new SpacingBetweenLines { Before = "60", After = "0", Line = "276", LineRule = LineSpacingRuleValues.Auto }),
                    new Run(new Text(v.GetString() ?? ""))));
            }
        }
    }

    // ── 段落 ──────────────────────────────────────────────────────
    private static void AddHeading(Body body, JsonElement blk)
    {
        var level = blk.TryGetProperty("level", out var l) ? l.GetInt32() : 1;
        var text = blk.GetProperty("text").GetString() ?? "";
        body.AppendChild(new Paragraph(
            new ParagraphProperties(new ParagraphStyleId { Val = $"Heading{level}" }),
            new Run(new Text(text))));
    }

    private static void AddParagraph(Body body, JsonElement blk)
    {
        var p = new Paragraph(new ParagraphProperties(new ParagraphStyleId { Val = "Normal" }));
        // 可选的行内富文本：runs 数组 [{text, bold?, italic?}]
        if (blk.TryGetProperty("runs", out var runs))
        {
            foreach (var r in runs.EnumerateArray())
            {
                var txt = r.GetProperty("text").GetString() ?? "";
                var run = new Run();
                if (r.TryGetProperty("bold", out var b) && b.GetBoolean())
                    run.Append(new RunProperties(new Bold()));
                if (r.TryGetProperty("italic", out var it) && it.GetBoolean())
                    run.Append(new RunProperties(new Italic()));
                run.Append(new Text(txt) { Space = SpaceProcessingModeValues.Preserve });
                p.Append(run);
            }
        }
        else
        {
            p.Append(new Run(new Text(blk.GetProperty("text").GetString() ?? "")
            { Space = SpaceProcessingModeValues.Preserve }));
        }
        body.AppendChild(p);
    }

    private static void AddBullet(Body body, JsonElement blk)
    {
        foreach (var item in blk.GetProperty("items").EnumerateArray())
        {
            body.AppendChild(new Paragraph(
                new ParagraphProperties(
                    new ParagraphStyleId { Val = "Normal" },
                    new Indentation { Left = "425", Hanging = "283" },
                    new SpacingBetweenLines { After = "0", Line = LineTwips.ToString(), LineRule = LineSpacingRuleValues.Auto }),
                new Run(new Text("• " + (item.GetString() ?? "")))));
        }
    }

    // ── 表格 ──────────────────────────────────────────────────────
    private static void AddTable(Body body, JsonElement blk, string baseDir)
    {
        if (blk.TryGetProperty("caption", out var cap))
        {
            body.AppendChild(new Paragraph(
                new ParagraphProperties(new ParagraphStyleId { Val = "Caption" },
                                        new KeepNext(),
                                        new Justification { Val = JustificationValues.Left }),
                new Run(new RunProperties(new Bold()), new Text(cap.GetString() ?? ""))));
        }

        var headers = blk.GetProperty("headers").EnumerateArray().Select(x => x.GetString() ?? "").ToArray();
        var rows = blk.GetProperty("rows").EnumerateArray()
            .Select(r => r.EnumerateArray().Select(x => x.GetString() ?? "").ToArray()).ToArray();
        var nCols = headers.Length;

        var tbl = new Table();
        var tblPr = new TableProperties();
        tblPr.Append(new TableWidth { Width = "5000", Type = TableWidthUnitValues.Pct });
        // 三线表（booktabs）：只有顶线、表头下横线、底线三条。
        // 无竖线、无行间横线、无底纹——行靠留白分隔，不靠框线。
        // 边框单位为 1/8 pt：Size=16 → 2 pt（顶/底线），Size=8 → 1 pt（表头线）。
        tblPr.Append(new TableBorders(
            new TopBorder { Val = BorderValues.Single, Size = 16U, Space = 0U },
            new BottomBorder { Val = BorderValues.Single, Size = 16U, Space = 0U },
            new LeftBorder { Val = BorderValues.None },
            new RightBorder { Val = BorderValues.None },
            new InsideHorizontalBorder { Val = BorderValues.None },
            new InsideVerticalBorder { Val = BorderValues.None }));
        tblPr.Append(new TableCellMarginDefault(
            new TopMargin { Width = "40", Type = TableWidthUnitValues.Dxa },
            new BottomMargin { Width = "40", Type = TableWidthUnitValues.Dxa },
            new TableCellLeftMargin { Width = 80, Type = TableWidthValues.Dxa },
            new TableCellRightMargin { Width = 80, Type = TableWidthValues.Dxa }));
        tbl.Append(tblPr);

        // 网格：正文宽度 11906 - 2*1440 = 9026 DXA
        const int totalW = 9026;
        // 可选的列宽比例（blk.widths，如 [2.0, 1, 1, ...]）；缺省等分。
        // 首列常是长标签列，等分会导致逐词折行，必须显式加宽。
        var ratios = new double[nCols];
        var hasW = blk.TryGetProperty("widths", out var wArr);
        if (hasW)
        {
            var wj = wArr.EnumerateArray().ToArray();
            for (var c = 0; c < nCols; c++)
                ratios[c] = c < wj.Length ? wj[c].GetDouble() : 1.0;
        }
        else
        {
            for (var c = 0; c < nCols; c++) ratios[c] = 1.0;
        }
        var rsum = ratios.Sum();
        var colW = new int[nCols];
        var acc = 0;
        for (var c = 0; c < nCols; c++)
        {
            colW[c] = c == nCols - 1 ? totalW - acc : (int)(totalW * ratios[c] / rsum);
            acc += colW[c];
        }
        var grid = new TableGrid();
        for (var c = 0; c < nCols; c++) grid.Append(new GridColumn { Width = colW[c].ToString() });
        tbl.Append(grid);

        // 表头行（跨页重复）。表头下的 1 pt 横线是三线表的中间那条，
        // 只能逐单元格设——表级 InsideHorizontal 已按三线表要求关闭。
        var hdr = new TableRow(new TableRowProperties(new TableHeader(), new CantSplit()));
        for (var c = 0; c < nCols; c++)
            hdr.Append(MakeCell(headers[c], colW[c], bold: true, ruleUnder: true));
        tbl.Append(hdr);

        foreach (var r in rows)
        {
            var tr = new TableRow(new TableRowProperties(new CantSplit()));
            for (var c = 0; c < nCols; c++)
                tr.Append(MakeCell(c < r.Length ? r[c] : "", colW[c], bold: false, ruleUnder: false));
            tbl.Append(tr);
        }
        body.AppendChild(tbl);
        body.AppendChild(new Paragraph(new ParagraphProperties(
            new ParagraphStyleId { Val = "Normal" },
            new SpacingBetweenLines { Before = "0", After = "0", Line = "240", LineRule = LineSpacingRuleValues.Auto })));
    }

    private static TableCell MakeCell(string text, int widthDxa, bool bold, bool ruleUnder)
    {
        var tcPr = new TableCellProperties(
            new TableCellWidth { Width = widthDxa.ToString(), Type = TableWidthUnitValues.Dxa });
        if (ruleUnder)
            tcPr.Append(new TableCellBorders(
                new BottomBorder { Val = BorderValues.Single, Size = 8U, Space = 0U }));
        // 单元格内必须至少有一个段落，否则 Word 判定文件损坏
        var para = new Paragraph(new ParagraphProperties(
            new ParagraphStyleId { Val = "TableText" },
            new Justification { Val = JustificationValues.Left },
            new SpacingBetweenLines { Before = "40", After = "40", Line = "240", LineRule = LineSpacingRuleValues.Auto }));
        var run = new Run();
        if (bold) run.Append(new RunProperties(new Bold()));
        run.Append(new Text(text) { Space = SpaceProcessingModeValues.Preserve });
        para.Append(run);
        return new TableCell(tcPr, para);
    }

    // ── 图 ────────────────────────────────────────────────────────
    private static void AddFigure(Body body, JsonElement blk, string baseDir)
    {
        if (blk.TryGetProperty("path", out var p))
        {
            var rel = Path.GetFullPath(Path.Combine(baseDir, p.GetString() ?? ""));
            if (File.Exists(rel))
            {
                var main = GetMainPart(body);
                var ip = main.AddImagePart(ImagePartType.Png);
                using (var fs = File.OpenRead(rel)) ip.FeedData(fs);
                var relId = main.GetIdOfPart(ip);
                var mm = blk.TryGetProperty("width_mm", out var w) ? w.GetDouble() : 160.0;
                var cx = (int)(mm * 36000);       // 1 mm = 36000 EMU
                var cy = (int)(cx * 0.62);
                var dpId = (UInt32)(++_drawingId);
                body.AppendChild(new Paragraph(
                    new ParagraphProperties(
                        new ParagraphStyleId { Val = "Normal" },
                        new Justification { Val = JustificationValues.Center },
                        new KeepNext(),
                        new SpacingBetweenLines { Before = "120", After = "60", Line = "240", LineRule = LineSpacingRuleValues.Auto }),
                    new Run(new Drawing(
                        new DocumentFormat.OpenXml.Drawing.Wordprocessing.Inline(
                            new DocumentFormat.OpenXml.Drawing.Wordprocessing.Extent { Cx = cx, Cy = cy },
                            new DocumentFormat.OpenXml.Drawing.Wordprocessing.EffectExtent { LeftEdge = 0L, TopEdge = 0L, RightEdge = 0L, BottomEdge = 0L },
                            new DocumentFormat.OpenXml.Drawing.Wordprocessing.DocProperties { Id = dpId, Name = Path.GetFileName(rel) },
                            new DocumentFormat.OpenXml.Drawing.Graphic(
                                new DocumentFormat.OpenXml.Drawing.GraphicData(
                                    new DocumentFormat.OpenXml.Drawing.Pictures.Picture(
                                        new DocumentFormat.OpenXml.Drawing.Pictures.NonVisualPictureProperties(
                                            new DocumentFormat.OpenXml.Drawing.Pictures.NonVisualDrawingProperties { Id = 0U, Name = Path.GetFileName(rel) },
                                            new DocumentFormat.OpenXml.Drawing.Pictures.NonVisualPictureDrawingProperties()),
                                        new DocumentFormat.OpenXml.Drawing.Pictures.BlipFill(
                                            new DocumentFormat.OpenXml.Drawing.Blip { Embed = relId },
                                            new DocumentFormat.OpenXml.Drawing.Stretch(new DocumentFormat.OpenXml.Drawing.FillRectangle())),
                                        new DocumentFormat.OpenXml.Drawing.Pictures.ShapeProperties(
                                            new DocumentFormat.OpenXml.Drawing.Transform2D(
                                                new DocumentFormat.OpenXml.Drawing.Offset { X = 0L, Y = 0L },
                                                new DocumentFormat.OpenXml.Drawing.Extents { Cx = cx, Cy = cy }),
                                            new DocumentFormat.OpenXml.Drawing.PresetGeometry(
                                                new DocumentFormat.OpenXml.Drawing.AdjustValueList())
                                            { Preset = DocumentFormat.OpenXml.Drawing.ShapeTypeValues.Rectangle })))
                                { Uri = "http://schemas.openxmlformats.org/drawingml/2006/picture" }))
                            { DistanceFromTop = 0U, DistanceFromBottom = 0U, DistanceFromLeft = 0U, DistanceFromRight = 0U }))));
            }
            else
            {
                Console.Error.WriteLine($"warn: figure not found: {rel}");
            }
        }
        if (blk.TryGetProperty("caption", out var cap))
        {
            body.AppendChild(new Paragraph(
                new ParagraphProperties(new ParagraphStyleId { Val = "Caption" },
                                        new Justification { Val = JustificationValues.Left }),
                new Run(new Text(cap.GetString() ?? ""))));
        }
    }

    private static uint _drawingId = 1U;

    private static MainDocumentPart GetMainPart(Body body)
        => body.Ancestors<Document>().First().MainDocumentPart!;
}
