#pragma once

#include <map>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <sstream>

#include <boost/version.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace libmpdataxx
{
  namespace output
  {
    namespace detail
    {
      template<int dim>
      class xdmf_writer
      {
        using ptree = boost::property_tree::ptree;
#if (BOOST_VERSION >= 105600)
        using xml_writer_settings = boost::property_tree::xml_writer_settings<std::string>;
#else
        using xml_writer_settings = boost::property_tree::xml_writer_settings<char>;
#endif

        struct data_item
        {
          blitz::TinyVector<int, dim> dimensions;
          static const std::string number_type;
          static const std::string format;
          std::string data;
          void add(ptree& node)
          {
            std::stringstream ss;
            for (auto d : dimensions)
              ss << d << ' ';
            ptree& dat_node = node.add("DataItem", data);
            dat_node.put("<xmlattr>.Dimensions", ss.str());
            dat_node.put("<xmlattr>.NumberType", number_type);
            dat_node.put("<xmlattr>.Format", format);
          }
        };

        struct topology
        {
          static const std::string topology_type;
          blitz::TinyVector<int, dim> dimensions;
          void add(ptree& node)
          {
            node.put("Topology.<xmlattr>.TopologyType", topology_type);
            std::stringstream ss;
            for (auto d : dimensions)
              ss << d << ' ';
            node.put("Topology.<xmlattr>.Dimensions", ss.str());
          }
        };

        struct geometry
        {
          static const std::string geometry_type;
          std::array<data_item, dim> coords;
          void add(ptree& node)
          {
            ptree& geo_node = node.add("Geometry", "");
            geo_node.put("<xmlattr>.GeometryType", geometry_type);
            for (auto &c : coords)
              c.add(geo_node);
          }
        };

        struct attribute
        {
          std::string name;
          static const std::string attribute_type;
          static const std::string center;
          mutable data_item item;
          void add(ptree& node)
          {
            ptree& attr_node = node.add("Attribute", "");
            attr_node.put("<xmlattr>.Name", name);
            attr_node.put("<xmlattr>.AttributeType", attribute_type);
            attr_node.put("<xmlattr>.Center", center);
            item.add(attr_node);
          }

          // to allow storing attributes in std::set
          friend bool operator<(const attribute &lhs, const attribute &rhs)
          {
            return lhs.name < rhs.name;
          }
        };

        static const std::string name;
        static const std::string grid_type;
        topology top;
        geometry geo;
        std::set<attribute> attrs;
        std::set<attribute> c_attrs;

        attribute make_attribute(const std::string& name,
                                 const blitz::TinyVector<int, dim>& dimensions)
        {
          attribute a;
          a.name = name;
          a.item.dimensions = dimensions - 1;
          return a;
        }

        public:

        void setup(const std::string& hdf_name,
                   const std::map<int, std::string>& dim_names,
                   const std::vector<std::string>& attr_names,
                   const blitz::TinyVector<int, dim>& dimensions)
        {
          top.dimensions = dimensions;

          for (const auto& dn : dim_names)
          {
            geo.coords[dn.first].dimensions = dimensions;
            geo.coords[dn.first].data = hdf_name + ":/" + dn.second;
          }

          for (const auto& n : attr_names)
          {
            attrs.insert(make_attribute(n, dimensions));
          }
        }
        

        void add_attribute(const std::string& name,
                                 const std::string& hdf_name,
                                 const blitz::TinyVector<int, dim>& dimensions)
        {
          attribute a = make_attribute(name, dimensions);
          a.item.data = hdf_name + ":/" + a.name;
          attrs.insert(a);
        }

        void add_const_attribute(const std::string& name,
                                 const std::string& hdf_name,
                                 const blitz::TinyVector<int, dim>& dimensions)
        {
          attribute a = make_attribute(name, dimensions);
          a.item.data = hdf_name + ":/" + a.name;
          c_attrs.insert(a);
        }

        void write(const std::string& xmf_name, const std::string& hdf_name, const double time)
        {
          for (auto& a : attrs)
          {
            a.item.data = hdf_name + ":/" + a.name;
          }

          ptree pt;
          ptree& grid_node = pt.put("Xdmf.Domain.Grid", "");
          grid_node.put("<xmlattr>.Name", name);
          grid_node.put("<xmlattr>.GridType", grid_type);
          grid_node.put("<xmlattr>.xml:id", "gid");
          grid_node.put("Time.<xmlattr>.Value", std::to_string(time));

          top.add(grid_node);

          geo.add(grid_node);

          for (auto a : attrs)
            a.add(grid_node);
          
          for (auto ca : c_attrs)
            ca.add(grid_node);

          xml_writer_settings settings('\t', 1);
          write_xml(xmf_name, pt, std::locale(), settings);
        }

        void write_temporal(const std::string& xmf_name, const std::vector<std::string>& timesteps)
        {

          ptree pt;
          ptree& xdmf_node = pt.put("Xdmf", "");
          xdmf_node.put("<xmlattr>.xmlns:xi", "http://www.w3.org/2001/XInclude");
          ptree& grid_node = xdmf_node.put("Domain.Grid", "");
          grid_node.put("<xmlattr>.Name", "TimeGrid");
          grid_node.put("<xmlattr>.GridType", "Collection");
          grid_node.put("<xmlattr>.CollectionType", "Temporal");

          for (auto ts : timesteps)
          {
            ptree& ts_node = grid_node.add("xi:include", "");
            ts_node.put("<xmlattr>.href", ts);
            ts_node.put("<xmlattr>.xpointer", "gid");
          }

          xml_writer_settings settings('\t', 1);
          write_xml(xmf_name, pt, std::locale(), settings);
        }

      };

      // compile time constants
      template<int dim>
      const std::string xdmf_writer<dim>::data_item::number_type = "Float";
      template<int dim>
      const std::string xdmf_writer<dim>::data_item::format = "HDF";
      template<>
      const std::string xdmf_writer<3>::topology::topology_type = "3DSMesh";
      template<>
      const std::string xdmf_writer<2>::topology::topology_type = "2DSMesh";
      template<>
      const std::string xdmf_writer<3>::geometry::geometry_type = "X_Y_Z";
      template<>
      const std::string xdmf_writer<2>::geometry::geometry_type = "X_Y";
      template<int dim>
      const std::string xdmf_writer<dim>::attribute::attribute_type = "Scalar";
      template<int dim>
      const std::string xdmf_writer<dim>::attribute::center = "Cell";
      template<int dim>
      const std::string xdmf_writer<dim>::name = "Grid";
      template<int dim>
      const std::string xdmf_writer<dim>::grid_type = "Uniform";
    }
  }
}
