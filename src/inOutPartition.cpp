#include <iostream>
#include <fstream>
#include <algorithm>
#include "BSP.h"

#include "graph_cut/GCoptimization.h"
#include "graph_cut/LinkedBlockList.cpp"
#include "graph_cut/GCoptimization.cpp"

// Label = 1 means INTERNAL
// Label = 0 means EXTERNAL

// Take first coplanar constraint associated to this face
// and return TRUE if the cell vertices are 'below' such constraint.
// bool isFirstConnCellBelowFace(BSPface& f, BSPcomplex* cpx, bool& is_degenerate)
bool isFirstConnCellBelowFace(BSPface& f, BSPcomplex* cpx)
{
    // is_degenerate = false;
    const uint32_t* cid = cpx->constraints_vrts.data() + f.coplanar_constraints[0] * 3;
    // std::cout << "cpx->constraints_vrts.data(): " << cpx->constraints_vrts.data()[0] << ", " << cpx->constraints_vrts.data()[1] << ", " << cpx->constraints_vrts.data()[2] << std::endl;
    // std::cout << "f.coplanar_constraints[0]: " << f.coplanar_constraints[0] << std::endl;
    // std::cout << "cid: " << cid[0] << ", " << cid[1] << ", " << cid[2] << std::endl;
    const genericPoint* pv1 = cpx->vertices[cid[0]];
    const genericPoint* pv2 = cpx->vertices[cid[1]];
    const genericPoint* pv3 = cpx->vertices[cid[2]];
    // std::cout << "pv1: " << *pv1 << std::endl;
    // std::cout << "pv2: " << *pv2 << std::endl;
    // std::cout << "pv3: " << *pv3 << std::endl;

    BSPcell& cell = cpx->cells[f.conn_cells[0]];
    // std::cout << "f.conn_cells[0]: " << f.conn_cells[0] << std::endl;
    uint64_t num_cellEdges = UINT64_MAX;
    uint32_t num_cellVrts = cpx->count_cellVertices(cell, &num_cellEdges);
    // std::cout << "num_cellEdges: " << num_cellEdges << ", num_cellVrts: " << num_cellVrts << std::endl;
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    cpx->list_cellVertices(cell, num_cellEdges, cell_vrts);
    // std::cout << "0: "  << std::endl;
    // std::cout << "cell_vrts: " << cell_vrts.size() << std::endl;

    for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 1;
    // std::cout << "1: "  << std::endl;

    // for (uint32_t vi : cell_vrts) {
    //     // std::cout << "2: " << std::endl;
    //     // const genericPoint* cv = cpx->vertices[vi];
    //     // std::cout << "cpx->vertices[vi]: " << vi << std::endl;
    //     // std::cout << "cv: " << *cv << std::endl;
    //     // std::cout << "cpx->vertices[vi]: " << vi << ", cv: " << *cv << ", vrts_visit: " << cpx->vrts_visit[vi] << std::endl;
    // }

    // std::cout << "- - - "  << std::endl;
    for (uint32_t vi : cell_vrts) if (!cpx->vrts_visit[vi])
    {
        const genericPoint* cv = cpx->vertices[vi];
        // std::cout << "cv: " << *cv << std::endl;
        const int o = genericPoint::orient3D(*cv, *pv1, *pv2, *pv3);
        // std::cout << "o: " << o << std::endl;
        if (o)
        {
            for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 0;
            return (o > 0);
        } else {
            // std::cout << "o: " << o << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
        }
    }

    // is_degenerate = true;
    ip_error("Degenerate cell\n");
    return false;
}

bool isSecondConnCellBelowFace(BSPface& f, BSPcomplex* cpx, bool& is_degenerate)
{
    is_degenerate = false;
    const uint32_t* cid = cpx->constraints_vrts.data() + f.coplanar_constraints[0] * 3;
    const genericPoint* pv1 = cpx->vertices[cid[0]];
    const genericPoint* pv2 = cpx->vertices[cid[1]];
    const genericPoint* pv3 = cpx->vertices[cid[2]];

    BSPcell& cell = cpx->cells[f.conn_cells[1]];
    uint64_t num_cellEdges = UINT64_MAX;
    uint32_t num_cellVrts = cpx->count_cellVertices(cell, &num_cellEdges);
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    cpx->list_cellVertices(cell, num_cellEdges, cell_vrts);

    for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 1;

    for (uint32_t vi : cell_vrts) if (!cpx->vrts_visit[vi])
    {
        const genericPoint* cv = cpx->vertices[vi];
        const int o = genericPoint::orient3D(*cv, *pv1, *pv2, *pv3);
        if (o)
        {
            for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 0;
            return (o > 0);
        } else {
        }
    }

    is_degenerate = true;
    ip_error("Degenerate cell\n");
    return false;
}

#define DETERMINANT3X3(a11, a12, a13, a21, a22, a23, a31, a32, a33) ((a11)*((a22)*(a33) - (a23)*(a32)) - (a12)*((a21)*(a33) - (a23)*(a31)) + (a13)*((a21)*(a32) - (a22)*(a31)))

double approxFaceArea(BSPface& face, BSPcomplex* cpx, const std::vector<double>& approxCoords)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);

    const double* acp = approxCoords.data();

    const double* tv0, * tv1, * tv2;
    double a = 0.0;
    tv0 = acp + vs[0] * 3;
    for (size_t i = 2; i < vs.size(); i++)
    {
        tv1 = acp + vs[i - 1] * 3;
        tv2 = acp + vs[i    ] * 3;
        a += DETERMINANT3X3(tv0[0], tv0[1], tv0[2], tv1[0], tv1[1], tv1[2], tv2[0], tv2[1], tv2[2]);
    }

    return fabs(a);
}


// Returns TRUE if face is part of the skin according to skin_colour
inline bool isSkinFace(const BSPface& face, uint32_t skin_colour)
{
    return face.colour & skin_colour;
}

inline void setInternalCell(BSPcell& c, uint32_t internal_label)
{
    c.place |= internal_label;
}

inline void setExternalCell(BSPcell& c)
{
    c.place = EXTERNAL;
}


void BSPcomplex::markInternalCells(uint32_t skin_colour, uint32_t internal_label, const std::vector<double>& face_areas)
{
    // Allocate dual graph: num cells + 1 to account for the external "ghost" cell
    GCoptimizationGeneralGraph gc((GCoptimization::SiteID)cells.size() + 1, 2);

    // gc is the dual graph of the cell complex
    // - a node in gc corresponds to a cell in the complex
    // - an arc in gc exists if two cells share a WHITE face

    // In gc, cells that share a white face are connected by an arc weighted on face area
    // The 'data cost' associated to each cell is:
    //  - total area of BLACK faces 'consistently oriented' with the cell, if label is EXTERNAL
    //  - 0, if label is INTERNAL
    //
    // The 'smooth cost' associated to each WHITE face is:
    //  - the face area, if label_A and label_B are different
    //  - 0 otherwise
    //
    // Note: a BLACK face is considered to be consistently oriented with one of its incident cells
    // if the cell is 'below' the first of its coplanar constraints.

    // evs == 1 if edge is on boundary of skin
    std::vector<uint8_t> evs(edges.size(), 0);
    for (BSPface& f : faces) if (isSkinFace(f, skin_colour))
      for (uint64_t eid : f.edges) if (evs[eid] < 2) evs[eid]++;

    // vvs == 1 if vertex is on boundary of skin
    std::vector<uint8_t> vvs(vertices.size(), 0);
    for (size_t i = 0; i < edges.size(); i++) if (evs[i] == 1)
    {
        const BSPedge& e = edges[i];
        vvs[e.vertices[0]] = vvs[e.vertices[1]] = 1;
    }

    std::vector<double> cell_costs_external(cells.size() + 1, 0.0);
    std::vector<double> cell_costs_internal(cells.size() + 1, 0.0);

    // std::cout << "BSPcomplex::markInternalCells" << std::endl;
    // std::cout << "cells.size(): " << cells.size() << std::endl;
    // std::cout << "edges.size(): " << edges.size() << ", vertices.size(): " << vertices.size() << ", faces.size(): " << faces.size() << std::endl;
    int tmp_cnt_ =0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);

        // if (isSkinFace(f, skin_colour))
        if (isSkinFace(f, skin_colour) && !(f.colour==RED))
        {
            // std::cout << "- - - - - - --  start - -- - - -- - -- " << std::endl;
            // std::cout << "f.conn_cells[1]: " << f.conn_cells[1] << ", cells.size(): " << cells.size() << std::endl;
            // std::cout << "i: " << i << ", f.colour: " << f.colour << std::endl;
            // std::cout << "cell1: " << cell1 << ", cell2: " << cell2 << std::endl;

            bool is_first_conn_cell_below_face;
            if ((in_out_label[cell1]!=-1) && in_out_label[cell2]!=-1) {
                if ((in_out_label[cell1]==0) && (in_out_label[cell2]==1)) { // cell1: in, cell2: out
                    is_first_conn_cell_below_face = true;
                } else if ((in_out_label[cell1]==1) && (in_out_label[cell2]==0)) { // cell1: out, cell2: in
                    is_first_conn_cell_below_face = false;
                } else {
                    std::cout << "WRONG IN/OUT LABELING! cell1: " << cell1 << ", "<< in_out_label[cell1] <<" | cell2: " << cell2 << ", "<< in_out_label[cell2] << std::endl;
                    continue;
                    is_first_conn_cell_below_face = isFirstConnCellBelowFace(f, this);
                    std::cout << "calculated: " << is_first_conn_cell_below_face << std::endl;
                }
                // if (is_first_conn_cell_below_face!=isFirstConnCellBelowFace(f, this)){tmp_cnt_++;}
                // std::cout << "labeled : " << is_first_conn_cell_below_face << ", calculated: "<< isFirstConnCellBelowFace(f, this) << std::endl;
            } else {
                std::cout << "NO INITIALIZED IN/OUT LABEL cell1: " << cell1 << ", cell2: " << cell2 << std::endl;
                is_first_conn_cell_below_face = isFirstConnCellBelowFace(f, this);
            }


            // const uint32_t* cid = constraints_vrts.data() + f.coplanar_constraints[0] * 3;
            // if ((cid[2]==37) || (cid[2]==39) || (cid[2]==38)) {
            // // if ((cid[0]==38) || (cid[1]==38) || (cid[2]==38)) {
            //     std::cout << "cid: " << cid[0] << ", " << cid[1] << ", " << cid[2] << std::endl;
            //     std::cout << "face_areas: " << face_areas[i] << std::endl;
            // }
            if (is_first_conn_cell_below_face)
            {
                cell_costs_external[cell1] += face_areas[i];
                cell_costs_internal[cell2] += face_areas[i];
            }
            else
            {
                cell_costs_external[cell2] += face_areas[i];
                cell_costs_internal[cell1] += face_areas[i];
            }

            // bool is_degenerate;
            // bool is_first_conn_cell_below_face = isFirstConnCellBelowFace(f, this, is_degenerate);

            // if (!is_degenerate) {
            //     if (is_first_conn_cell_below_face)
            //     {
            //         cell_costs_external[cell1] += face_areas[i];
            //         cell_costs_internal[cell2] += face_areas[i];
            //     }
            //     else
            //     {
            //         cell_costs_external[cell2] += face_areas[i];
            //         cell_costs_internal[cell1] += face_areas[i];
            //     }
            // } else {
            //     bool is_second_conn_cell_below_face = isSecondConnCellBelowFace(f, this, is_degenerate);
            //     if (!is_degenerate) {
            //         if (is_first_conn_cell_below_face)
            //         {
            //             cell_costs_external[cell2] += face_areas[i];
            //             cell_costs_internal[cell1] += face_areas[i];
            //         }
            //         else
            //         {
            //             cell_costs_external[cell1] += face_areas[i];
            //             cell_costs_internal[cell2] += face_areas[i];
            //         }
            //     }
            // }

            // if (is_first_conn_cell_below_face)

            // if (isFirstConnCellBelowFace(f, this))
            // {
            //     cell_costs_external[cell1] += face_areas[i];
            //     cell_costs_internal[cell2] += face_areas[i];
            // }
            // else
            // {
            //     cell_costs_external[cell2] += face_areas[i];
            //     cell_costs_internal[cell1] += face_areas[i];
            // }
            // std::cout << "- - - - - - -- - end -- - - -- - -- " << std::endl;
        }
    }
    // std::cout << "mis labeld num: " <<tmp_cnt_ << std::endl;
    // std::cout << "set cell_costss done" << std::endl;
    int tmp_cnt = 0;
    int already_cnt = 0;
    int notyet_cnt = 0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        if (!isSkinFace(f, skin_colour) || (f.colour==RED))
        {
            // 'w' is an additional weight for arcs that promotes the cut of arcs corresp. to faces having
            // all their vertices on the boundary of the input surface (i.e. promote hole filling)
            double w = 0.1;
            // double w_extra = 10.0; // jywq test
            for (uint64_t eid : f.edges) if (vvs[edges[eid].vertices[0]] == 0 || vvs[edges[eid].vertices[1]] == 0) { w = 1.0; break; }
            if ((f.colour==RED)) {
                if(w==1.0) {
                    already_cnt++;
                } else {
                    notyet_cnt++;
                }
                w=10.0;
                tmp_cnt++;
            }
            const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);
            gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, face_areas[i]*w);
            // gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, face_areas[i]*w * w_extra);

            // 禁止更改cell都在内部的面
            if ((in_out_label[cell1]==0) && (in_out_label[cell2]==0)) { // cell1: in, cell2: in
                gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, 1.0);
            }
        }
    }
    // std::cout << "setNeighbors done" << std::endl;
    // std::cout << "!!!!!!!!!!! RED NUM: " << tmp_cnt << std::endl;
    // std::cout << "already: " << already_cnt<<std::endl;
    // std::cout << "notyet_cnt: " << notyet_cnt<<std::endl;

    const double int_weight = 0.1; // Internal cell penalization less than external to avoid artifacts at intersections
    // for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, cell_costs_external[i]);
    // for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, cell_costs_internal[i] * int_weight);
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, (cell_costs_external[i]>0) ? 1.0:0.0); // skin face must be selected
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, (cell_costs_internal[i]>0) ? 1.0:0.0); // skin face must be selected
    gc.setDataCost((GCoptimization::SiteID)cells.size(), 1, 1.0); // Ghost cell must be external
    // std::cout << "setDataCost done" << std::endl;

    // Run graph cut algorithm
    // I.e., label all the cells so that the total data cost + smooth cost is minimized
    gc.swap();
    // std::cout << "graph cut done" << std::endl;

    for (size_t i = 0; i < cells.size(); i++)
        if (gc.whatLabel((GCoptimization::SiteID)i)) setInternalCell(cells[i], internal_label);
}


void BSPcomplex::markInternalCells_distance(uint32_t skin_colour, uint32_t internal_label, const std::vector<double>& face_areas, const std::vector<double>& face_dists)
{
    // Allocate dual graph: num cells + 1 to account for the external "ghost" cell
    GCoptimizationGeneralGraph gc((GCoptimization::SiteID)cells.size() + 1, 2);

    // gc is the dual graph of the cell complex
    // - a node in gc corresponds to a cell in the complex
    // - an arc in gc exists if two cells share a WHITE face

    // In gc, cells that share a white face are connected by an arc weighted on face area
    // The 'data cost' associated to each cell is:
    //  - total area of BLACK faces 'consistently oriented' with the cell, if label is EXTERNAL
    //  - 0, if label is INTERNAL
    //
    // The 'smooth cost' associated to each WHITE face is:
    //  - the face area, if label_A and label_B are different
    //  - 0 otherwise
    //
    // Note: a BLACK face is considered to be consistently oriented with one of its incident cells
    // if the cell is 'below' the first of its coplanar constraints.

    // evs == 1 if edge is on boundary of skin
    std::vector<uint8_t> evs(edges.size(), 0);
    for (BSPface& f : faces) if (isSkinFace(f, skin_colour))
      for (uint64_t eid : f.edges) if (evs[eid] < 2) evs[eid]++;

    // vvs == 1 if vertex is on boundary of skin
    std::vector<uint8_t> vvs(vertices.size(), 0);
    for (size_t i = 0; i < edges.size(); i++) if (evs[i] == 1)
    {
        const BSPedge& e = edges[i];
        vvs[e.vertices[0]] = vvs[e.vertices[1]] = 1;
    }

    std::vector<double> cell_costs_external(cells.size() + 1, 0.0);
    std::vector<double> cell_costs_internal(cells.size() + 1, 0.0);

    // std::cout << "BSPcomplex::markInternalCells" << std::endl;
    // std::cout << "cells.size(): " << cells.size() << std::endl;
    // std::cout << "edges.size(): " << edges.size() << ", vertices.size(): " << vertices.size() << ", faces.size(): " << faces.size() << std::endl;
    int tmp_cnt_ =0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);

        if (isSkinFace(f, skin_colour)) // constraint faces
        {
            bool is_first_conn_cell_below_face;
            if ((in_out_label[cell1]!=-1) && in_out_label[cell2]!=-1) {
                if ((in_out_label[cell1]==0) && (in_out_label[cell2]==1)) { // cell1: in, cell2: out
                    is_first_conn_cell_below_face = true;
                } else if ((in_out_label[cell1]==1) && (in_out_label[cell2]==0)) { // cell1: out, cell2: in
                    is_first_conn_cell_below_face = false;
                } else {
                    std::cout << "WRONG IN/OUT LABELING! cell1: " << cell1 << ", "<< in_out_label[cell1] <<" | cell2: " << cell2 << ", "<< in_out_label[cell2] << std::endl;
                    continue;
                    is_first_conn_cell_below_face = isFirstConnCellBelowFace(f, this);
                    std::cout << "calculated: " << is_first_conn_cell_below_face << std::endl;
                }
            } else {
                std::cout << "NO INITIALIZED IN/OUT LABEL cell1: " << cell1 << ", cell2: " << cell2 << std::endl;
                is_first_conn_cell_below_face = isFirstConnCellBelowFace(f, this);
            }

            if (is_first_conn_cell_below_face)
            {
                // cell_costs_external[cell1] += face_areas[i];
                // cell_costs_internal[cell2] += face_areas[i];
                cell_costs_external[cell1] += std::abs(face_dists[i]);
                cell_costs_internal[cell2] += std::abs(face_dists[i]);
            }
            else
            {
                // cell_costs_external[cell2] += face_areas[i];
                // cell_costs_internal[cell1] += face_areas[i];
                cell_costs_external[cell2] += std::abs(face_dists[i]);
                cell_costs_internal[cell1] += std::abs(face_dists[i]);
            }
        }
    }
    // std::cout << "mis labeld num: " <<tmp_cnt_ << std::endl;
    // std::cout << "set cell_costss done" << std::endl;
    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        if (!isSkinFace(f, skin_colour)) // extra faces
        {
            // 'w' is an additional weight for arcs that promotes the cut of arcs corresp. to faces having
            // all their vertices on the boundary of the input surface (i.e. promote hole filling)
            double w = 0.1;
            for (uint64_t eid : f.edges) if (vvs[edges[eid].vertices[0]] == 0 || vvs[edges[eid].vertices[1]] == 0) { w = 1.0; break; }
            const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);
            gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, std::abs(face_dists[i]*w));
        }
    }
    // std::cout << "setNeighbors done" << std::endl;

    const double int_weight = 0.1; // Internal cell penalization less than external to avoid artifacts at intersections
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, cell_costs_external[i]);
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, cell_costs_internal[i] * int_weight);
    gc.setDataCost((GCoptimization::SiteID)cells.size(), 1, 1.0); // Ghost cell must be external
    // std::cout << "setDataCost done" << std::endl;

    // Run graph cut algorithm
    // I.e., label all the cells so that the total data cost + smooth cost is minimized
    gc.swap();
    // std::cout << "graph cut done" << std::endl;

    for (size_t i = 0; i < cells.size(); i++)
        if (gc.whatLabel((GCoptimization::SiteID)i)) setInternalCell(cells[i], internal_label);
}

#include <Eigen/Dense>
#include <cmath>

double calculateTriangleArea(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    // Calculate vectors representing two sides of the triangle
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;

    // Compute the cross product of the two vectors
    Eigen::Vector3d crossProduct = v1.cross(v2);

    // The area of the triangle is half the magnitude of the cross product
    double area = 0.5 * crossProduct.norm();

    return area;
}


double jywqFaceArea(BSPface& face, BSPcomplex* cpx, const std::vector<double>& approxCoords)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);

    Eigen::MatrixXd V(3,3);

    for (int j=0; j<3;j++) {
        V(j, 0) = approxCoords[vs[j] * 3];
        V(j, 1) = approxCoords[vs[j] * 3 + 1];
        V(j, 2) = approxCoords[vs[j] * 3 + 2];
    }
    double area = calculateTriangleArea(V.row(0), V.row(1), V.row(2));

    return area;
}

double jywqFaceDist(BSPface& face, BSPcomplex* cpx)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);

    double face_dist = 0;
    for (int j=0; j<3;j++) {
        face_dist += cpx->vert_dists[vs[j]];
    }

    return face_dist;
}

void BSPcomplex::constraintsSurface_complexPartition(bool two_files, bool distance_policy)
{
    // Make all cells external
    // std::cout << "distance_policy: "<< distance_policy << std::endl;
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);

    // Clear vrts_visit for use in isFirstConnCellBelowFace()
    for (size_t i = 0; i < vertices.size(); i++) vrts_visit[i] = 0;

    // Precalculate approximate vertex coordinates for use in approxFaceArea()
    std::vector<double> approxCoords(vertices.size() * 3);
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i * 3], approxCoords[i * 3 + 1], approxCoords[i * 3 + 2]);

    // Precalculate approximate face areas for use in markInternalCells()
    std::vector<double> face_areas(faces.size(), 0.0);
    std::vector<double> face_dists(faces.size(), 0.0);
    double tot_face_area = 0.0;
    double tot_face_dist = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        // face_areas[i] = approxFaceArea(faces[i], this, approxCoords);
        face_areas[i] = jywqFaceArea(faces[i], this, approxCoords);
        tot_face_area += std::abs(face_areas[i]);
        if ( distance_policy ) {
            // face_dists[i] = jywqFaceDist(faces[i], this);
            // tot_face_dist += std::abs(face_dists[i]);
            face_dists[i] = std::abs(jywqFaceDist(faces[i], this));
            tot_face_dist += face_dists[i];
        }

        // BSPface& f = faces[i];
        // // if (isSkinFace(f, skin_colour))
        // if (isSkinFace(f, BLACK_A) && !(f.colour==RED))
        // {
        //     const uint32_t* cid = constraints_vrts.data() + f.coplanar_constraints[0] * 3;
        //     if ((cid[2]==37) || (cid[2]==39) || (cid[2]==38)) {
        //     // if ((cid[0]==38) || (cid[1]==38) || (cid[2]==38)) {
        //         std::cout << "i: "<< i <<", cid: " << cid[0] << ", " << cid[1] << ", " << cid[2] << std::endl;
        //         std::vector<uint32_t> vs(f.edges.size(), 0);
        //         this->list_faceVertices(f, vs);
        //         const uint32_t* mv = f.meshVertices; 
        //         Eigen::MatrixXd V(3,3);
   
        //         for (int j=0; j<3;j++) {
        //             const genericPoint* gp = vertices[mv[j]];
        //             const explicitPoint3D& ep = gp->toExplicit3D();
        //             std::cout << "approxCoords " << vs[j]*3 << ": "<< approxCoords[vs[j]*3] << ", " << approxCoords[vs[j] * 3 + 1] << ", " << approxCoords[vs[j] * 3 + 2]<< std::endl;
        //             std::cout << "ori          " << mv[j] << ": "<< ep.X() << ", " << ep.Y() << ", " << ep.Z() << std::endl;
        //             V(j, 0) = ep.X();
        //             V(j, 1) = ep.Y();
        //             V(j, 2) = ep.Z();

        //         }
        //         double area = calculateTriangleArea(V.row(0), V.row(1), V.row(2));
        //         std::cout << "face_areas1: " << face_areas[i] << std::endl;
        //         std::cout << "face_areas2: " << area << std::endl;
        //         std::cout << "face_areas3: " << approxFaceArea(faces[i], this, approxCoords) << std::endl;
        //     }
        // }
    }

    // Normalize all areas to avoid overflows in graphcut
    for (size_t i = 0; i < faces.size(); i++) face_areas[i] /= tot_face_area;
    if ( distance_policy ) { for (size_t i = 0; i < faces.size(); i++) face_dists[i] /= tot_face_dist; }
    // std::cout << "tot_face_area: " << tot_face_area << std::endl;
    // std::cout << "tot_face_dist: " << tot_face_dist << std::endl;
    // std::cout << "distance_policy: "<< distance_policy<<", "<< (distance_policy ? "0":"1") << std::endl;

    // std::cout << "markInternalCells start" << std::endl;
    if (two_files)
    {
        markInternalCells(BLACK_A, INTERNAL_A, face_areas);
        markInternalCells(BLACK_B, INTERNAL_B, face_areas);
    }
    else markInternalCells(BLACK_A, INTERNAL_A, distance_policy ? face_dists:face_areas);
    // else markInternalCells(BLACK_A, INTERNAL_A, face_dists);
    // else distance_policy ? markInternalCells_distance(BLACK_A, INTERNAL_A, face_areas, face_dists) : markInternalCells(BLACK_A, INTERNAL_A, face_areas);
}

